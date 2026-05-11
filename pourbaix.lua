-- This file is part of the pourbaix package.
-- Version: 1.0 (2026-05-10).
-- License: LaTeX Project Public License v1.3c or later.
-- Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>.

-- pourbaix.lua
-- LuaLaTeX/pgfplots engine for potential-pH diagrams.
-- Core computational engine; corrosion-specific features are intentionally excluded.

local M = {}
M.data = require("pourbaix-data")

local function log10(x)
  return math.log(x) / math.log(10)
end

local function finite(x)
  return x == x and x ~= math.huge and x ~= -math.huge
end

local function number(x, default)
  local n = tonumber(x)
  if n == nil then return default end
  return n
end

local function trim(s)
  return tostring(s or ""):gsub("^%s+", ""):gsub("%s+$", "")
end

local function truthy(v, default)
  if v == nil or v == "" then return default end
  v = tostring(v):lower()
  if v == "true" or v == "yes" or v == "on" or v == "1" or v == "oui" then return true end
  if v == "false" or v == "no" or v == "off" or v == "0" or v == "non" then return false end
  return default
end

local function fmt(x)
  if not finite(x) then return tostring(x) end
  if math.abs(x) < 1e-12 then x = 0 end
  return string.format("%.10g", x)
end

local function fmt_sig(x, sig)
  if not finite(x) then return tostring(x) end
  sig = math.floor(number(sig, 6) or 6)
  if sig < 1 then sig = 1 end
  if sig > 15 then sig = 15 end
  if math.abs(x) < 1e-12 then x = 0 end
  return string.format("%." .. sig .. "g", x)
end

local function letter_name(n)
  n = math.floor(number(n, 1) or 1)
  if n < 1 then n = 1 end
  local s = ""
  while n > 0 do
    n = n - 1
    s = string.char(97 + (n % 26)) .. s
    n = math.floor(n / 26)
  end
  return s
end

local function upper_letter_name(n)
  return string.upper(letter_name(n))
end

local function csv(list)
  return table.concat(list, ",")
end

local function shallow_copy(t)
  local r = {}
  for k, v in pairs(t) do r[k] = v end
  return r
end

local function list_species(elementData)
  local ids = {}
  for id, _ in pairs(elementData.species) do ids[#ids+1] = id end
  table.sort(ids, function(a,b)
    local oa = elementData.species[a].oxState or 0
    local ob = elementData.species[b].oxState or 0
    if oa == ob then return a < b end
    return oa < ob
  end)
  return ids
end

function M.active_set(elementData, speciesList)
  local set = {}
  if speciesList == nil or speciesList == "" or speciesList == "all" then
    for id, _ in pairs(elementData.species) do set[id] = true end
    return set
  end
  for item in tostring(speciesList):gmatch("([^,]+)") do
    local id = trim(item)
    if id ~= "" and elementData.species[id] then set[id] = true end
  end
  return set
end

function M.normalizeConcentrationConvention(v)
  v = trim(v or "species"):lower()
  v = v:gsub("[_-]", " "):gsub("%s+", " ")
  if v == "" or v == "species" or v == "all species" or v == "solution species"
     or v == "species concentration"
     or v == "concentration species" then
    return "species"
  end
  if v == "element species equal" or v == "element species equality"
     or v == "species equal" or v == "equal species"

     or v == "element/species" or v == "total element species equal" then
    return "element_species_equal"
  end
  if v == "element element equal" or v == "element equality"
     or v == "elements equal" or v == "equal elements"
     or v == "element concentrations equal" or v == "element/element"
     or v == "total element element equal" then
    return "element_element_equal"
  end
  -- Unknown value: keep the historical behavior rather than stopping compilation.
  return "species"
end

local function metal_count(species)
  if not species then return 1 end
  local n = number(species.nMetal, 1)
  if not n or n <= 0 then return 1 end
  return n
end

function M.logActivity(species, cTrace, convention, partnerSpecies)
  if not species or species.phase ~= "aq" then return 0 end

  convention = M.normalizeConcentrationConvention(convention)
  if convention == "species" then
    return log10(cTrace)
  end

  local nu = metal_count(species)
  local cSpecies = cTrace / nu

  if partnerSpecies and partnerSpecies.phase == "aq" then
    local nuPartner = metal_count(partnerSpecies)
    if convention == "element_species_equal" then
      -- [A] = [B] = c_tot/(nu_A + nu_B).
      cSpecies = cTrace / (nu + nuPartner)
    elseif convention == "element_element_equal" then
      -- nu_A[A] = nu_B[B] = c_tot/2.
      cSpecies = cTrace / (2 * nu)
    end
  end

  if cSpecies <= 0 then cSpecies = cTrace end
  return log10(cSpecies)
end

function M.reduceAcidBaseChain(elementData, oxState, cTrace, activeSpeciesSet, convention)
  convention = M.normalizeConcentrationConvention(convention)
  local filtered = {}
  for _, eq in ipairs(elementData.acidBaseEquilibria or {}) do
    local sp1 = elementData.species[eq.sideAcid]
    local sp2 = elementData.species[eq.sideBase]
    if sp1 and sp2 and sp1.oxState == oxState and sp2.oxState == oxState
       and activeSpeciesSet[eq.sideAcid] and activeSpeciesSet[eq.sideBase] then
      filtered[#filtered+1] = eq
    end
  end
  if #filtered == 0 then return {} end

  local function computePHb(eq)
    local nA = eq.nAcid or 1
    local nB = eq.nBase or 1
    local baseSp = elementData.species[eq.sideBase]
    local acidSp = elementData.species[eq.sideAcid]
    local aBase = M.logActivity(baseSp, cTrace, convention, acidSp)
    local aAcid = M.logActivity(acidSp, cTrace, convention, baseSp)
    return (eq.pKa + nB * aBase - nA * aAcid) / eq.nH
  end

  local chain = {}
  for _, eq in ipairs(filtered) do
    local e = shallow_copy(eq)
    e.pHb = computePHb(e)
    chain[#chain+1] = e
  end

  -- User-provided equilibria may arrive in arbitrary order. Sort transitions
  -- by increasing boundary pH before reducing the chain; otherwise a chain
  -- such as H2S/HS-/S2- entered as HS-/S2-, H2S/HS- would be merged
  -- incorrectly into H2S/S2-.
  table.sort(chain, function(a,b)
    if math.abs((a.pHb or 0) - (b.pHb or 0)) < 1e-12 then
      return tostring(a.sideAcid or "") < tostring(b.sideAcid or "")
    end
    return (a.pHb or 0) < (b.pHb or 0)
  end)

  local changed = true
  while changed do
    changed = false
    for i = 1, #chain - 1 do
      if chain[i+1].pHb <= chain[i].pHb then
        local merged = {
          sideAcid = chain[i].sideAcid,
          sideBase = chain[i+1].sideBase,
          nAcid = 1,
          nBase = 1,
          nH = chain[i].nH + chain[i+1].nH,
          nH2O = (chain[i].nH2O or 0) + (chain[i+1].nH2O or 0),
          pKa = chain[i].pKa + chain[i+1].pKa,
        }
        merged.pHb = computePHb(merged)
        table.remove(chain, i+1)
        chain[i] = merged
        changed = true
        break
      end
    end
  end
  return chain
end

function M.dominantFormAtOxState(elementData, oxState, pH, cTrace, activeSpeciesSet, convention)
  local candidates = {}
  for id, sp in pairs(elementData.species) do
    if sp.oxState == oxState and activeSpeciesSet[id] then
      candidates[#candidates+1] = id
    end
  end
  table.sort(candidates)
  if #candidates == 0 then return nil end
  if #candidates == 1 then return candidates[1] end

  local boundaries = M.reduceAcidBaseChain(elementData, oxState, cTrace, activeSpeciesSet, convention)
  if #boundaries == 0 then return candidates[1] end

  if pH < boundaries[1].pHb then return boundaries[1].sideAcid end
  for i = 1, #boundaries - 1 do
    if pH >= boundaries[i].pHb and pH < boundaries[i+1].pHb then
      return boundaries[i].sideBase
    end
  end
  return boundaries[#boundaries].sideBase
end

function M.Eeq(couple, oxSpecies, redSpecies, pH, cTrace, convention)
  local nOx  = couple.nOx  or 1
  local nRed = couple.nRed or 1
  local nHR  = couple.nH_right or 0
  local effNH = (couple.nH or 0) - nHR
  local logAOx  = M.logActivity(oxSpecies, cTrace, convention, redSpecies)
  local logARed = M.logActivity(redSpecies, cTrace, convention, oxSpecies)
  return couple.E0
    - 0.06 * effNH / couple.n * pH
    + (0.06 / couple.n) * (nOx * logAOx - nRed * logARed)
end

function M.redoxLineCoefficients(couple, oxSpecies, redSpecies, cTrace, convention)
  local nOx  = couple.nOx  or 1
  local nRed = couple.nRed or 1
  local nHR  = couple.nH_right or 0
  local effNH = (couple.nH or 0) - nHR
  local logAOx  = M.logActivity(oxSpecies, cTrace, convention, redSpecies)
  local logARed = M.logActivity(redSpecies, cTrace, convention, oxSpecies)
  local slope = -0.06 * effNH / couple.n
  local intercept = couple.E0 + (0.06 / couple.n) * (nOx * logAOx - nRed * logARed)
  return slope, intercept
end

local function has_intermediate_oxidation_state(elementData, oxId, redId)
  local oxSp = elementData.species and elementData.species[oxId]
  local redSp = elementData.species and elementData.species[redId]
  if not oxSp or not redSp then return false end
  local hi, lo = oxSp.oxState, redSp.oxState
  if hi == nil or lo == nil then return false end
  if hi < lo then hi, lo = lo, hi end
  for id, sp in pairs(elementData.species or {}) do
    if id ~= oxId and id ~= redId and sp.oxState and sp.oxState > lo and sp.oxState < hi then
      return true
    end
  end
  return false
end

function M.findCouple(elementData, oxId, redId)
  for _, c in ipairs(elementData.redoxCouples or {}) do
    if c.ox == oxId and c.red == redId then
      -- Thermodynamic-cycle couples generated automatically between
      -- non-adjacent oxidation states document cycles, but they must not
      -- short-circuit intermediate oxidation states while the diagram is built.
      -- If the intermediate state is truly unstable, buildLadder reconstructs
      -- a virtual disproportionation boundary through _deriveImplicitEeqCouple.
      if c._autoE0 and has_intermediate_oxidation_state(elementData, oxId, redId) then
        -- skip
      else
        return c
      end
    end
  end
  return nil
end

function M._deriveImplicitEeqCouple(elementData, highId, prevId, ladder, pH, cTrace, convention, activeSpeciesSet)
  local highSp = elementData.species[highId]
  local prevSp = elementData.species[prevId]
  if not highSp or not prevSp then return nil end
  local n_hp = highSp.oxState - prevSp.oxState
  if n_hp <= 0 then return nil end

  local candidates, seen = {}, {}
  local function add_candidate(id)
    if id and id ~= highId and id ~= prevId and not seen[id] then
      if (not activeSpeciesSet) or activeSpeciesSet[id] then
        seen[id] = true
        candidates[#candidates+1] = id
      end
    end
  end
  for _, item in ipairs(ladder or {}) do add_candidate(item.form) end
  for id, _ in pairs(elementData.species or {}) do add_candidate(id) end

  local function build_virtual_via_lower(yForm, ySp, c1, c2)
    local nOxC1, nRedC1 = c1.nOx or 1, c1.nRed or 1
    local nOxC2, nRedC2 = c2.nOx or 1, c2.nRed or 1
    if not (nOxC1 == 1 and nRedC1 == 1 and nOxC2 == 1 and nRedC2 == 1) then return nil end
    local n_hy = highSp.oxState - ySp.oxState
    local n_py = prevSp.oxState - ySp.oxState
    local nH_L = (c1.nH or 0) + (c2.nH_right or 0)
    local nH_R = (c1.nH_right or 0) + (c2.nH or 0)
    local nW_L = (c1.nH2O_left or 0) + (c2.nH2O or 0)
    local nW_R = (c1.nH2O or 0) + (c2.nH2O_left or 0)
    local dnH = nH_L - nH_R
    local dnW = nW_L - nW_R
    local E0_hp = (n_hy * c1.E0 - n_py * c2.E0) / n_hp
    return {
      ox = highId, red = prevId, n = n_hp, E0 = E0_hp,
      nH = math.max(0, dnH), nH_right = math.max(0, -dnH),
      nH2O_left = math.max(0, dnW), nH2O = math.max(0, -dnW),
      _virtual = true, _viaY = yForm,
    }
  end

  local function build_virtual_via_between(yForm, ySp, c1, c2)
    local nOxC1, nRedC1 = c1.nOx or 1, c1.nRed or 1
    local nOxC2, nRedC2 = c2.nOx or 1, c2.nRed or 1
    if not (nOxC1 == 1 and nRedC1 == 1 and nOxC2 == 1 and nRedC2 == 1) then return nil end
    local n_hy = highSp.oxState - ySp.oxState
    local n_yp = ySp.oxState - prevSp.oxState
    local nH_L = (c1.nH or 0) + (c2.nH or 0)
    local nH_R = (c1.nH_right or 0) + (c2.nH_right or 0)
    local nW_L = (c1.nH2O_left or 0) + (c2.nH2O_left or 0)
    local nW_R = (c1.nH2O or 0) + (c2.nH2O or 0)
    local dnH = nH_L - nH_R
    local dnW = nW_L - nW_R
    local E0_hp = (n_hy * c1.E0 + n_yp * c2.E0) / n_hp
    return {
      ox = highId, red = prevId, n = n_hp, E0 = E0_hp,
      nH = math.max(0, dnH), nH_right = math.max(0, -dnH),
      nH2O_left = math.max(0, dnW), nH2O = math.max(0, -dnW),
      _virtual = true, _viaY = yForm,
    }
  end

  -- Historical case: the intermediate lies below the current reductant.
  -- Derive high/prev by subtracting prev/y from high/y.
  for _, yForm in ipairs(candidates) do
    local ySp = elementData.species[yForm]
    if ySp and ySp.oxState and ySp.oxState < prevSp.oxState then
      local c1 = M.findCouple(elementData, highId, yForm)
      local c2 = M.findCouple(elementData, prevId, yForm)
      if c1 and c2 then
        local vc = build_virtual_via_lower(yForm, ySp, c1, c2)
        if vc then return { Eeq = M.Eeq(vc, highSp, prevSp, pH, cTrace, convention), virtualCouple = vc } end
      end
    end
  end

  -- Additional case: the intermediate lies between oxidant and reductant.
  -- It may have been popped from the redox convex envelope; the direct
  -- high/prev boundary must then be reconstructed by adding high/y and y/prev.
  -- This is the sulfate/HS- then sulfate/S2- case in the sulfur diagram.
  for _, yForm in ipairs(candidates) do
    local ySp = elementData.species[yForm]
    if ySp and ySp.oxState and ySp.oxState > prevSp.oxState and ySp.oxState < highSp.oxState then
      local c1 = M.findCouple(elementData, highId, yForm)
      local c2 = M.findCouple(elementData, yForm, prevId)
      if c1 and c2 then
        local vc = build_virtual_via_between(yForm, ySp, c1, c2)
        if vc then return { Eeq = M.Eeq(vc, highSp, prevSp, pH, cTrace, convention), virtualCouple = vc } end
      end
    end
  end

  return nil
end

function M.buildLadder(elementData, pH, cTrace, activeSpeciesSet, convention)
  convention = M.normalizeConcentrationConvention(convention)
  local present = {}
  for id, _ in pairs(activeSpeciesSet) do
    local sp = elementData.species[id]
    if sp then present[sp.oxState] = true end
  end
  local oxStates = {}
  for n, _ in pairs(present) do oxStates[#oxStates+1] = n end
  table.sort(oxStates)
  if #oxStates == 0 then return {} end

  local formByOx = {}
  for _, n in ipairs(oxStates) do
    formByOx[n] = M.dominantFormAtOxState(elementData, n, pH, cTrace, activeSpeciesSet, convention)
  end

  local ladder = {{ form = formByOx[oxStates[1]], E_to_prev = -math.huge, couple = nil }}
  for i = 2, #oxStates do
    local high = formByOx[oxStates[i]]
    if high then
      while #ladder > 0 do
        local prev = ladder[#ladder]
        local couple = M.findCouple(elementData, high, prev.form)
        local Eeq = nil
        if couple then
          Eeq = M.Eeq(couple, elementData.species[high], elementData.species[prev.form], pH, cTrace, convention)
        else
          local derived = M._deriveImplicitEeqCouple(elementData, high, prev.form, ladder, pH, cTrace, convention, activeSpeciesSet)
          if derived then
            couple = derived.virtualCouple
            Eeq = derived.Eeq
          end
        end
        if couple and Eeq > prev.E_to_prev then
          ladder[#ladder+1] = { form = high, E_to_prev = Eeq, couple = couple }
          break
        end
        if #ladder == 1 then break end
        ladder[#ladder] = nil
      end
    end
  end
  return ladder
end

function M.dominantSpeciesAt(elementData, pH, E, cTrace, activeSpeciesSet, convention)
  local ladder = M.buildLadder(elementData, pH, cTrace, activeSpeciesSet, convention)
  if #ladder == 0 then return nil end
  local currentForm = ladder[1].form
  for i = 2, #ladder do
    if E > ladder[i].E_to_prev then currentForm = ladder[i].form else break end
  end
  return currentForm
end

function M.default_range(elementData, range)
  range = range or {}
  local er = elementData.eRange or {-1.0, 2.0}
  return {
    pHmin = number(range.pHmin, 0),
    pHmax = number(range.pHmax, 14),
    Emin = number(range.Emin, er[1]),
    Emax = number(range.Emax, er[2]),
  }
end

local function split_point_groups(points, dpH)
  local groups = {}
  local current = {}
  local prevx = nil
  for _, p in ipairs(points) do
    if prevx and math.abs(p[1] - prevx) > 1.75 * dpH then
      if #current > 1 then groups[#groups+1] = current end
      current = {}
    end
    current[#current+1] = p
    prevx = p[1]
  end
  if #current > 1 then groups[#groups+1] = current end
  return groups
end

local function couple_signature(couple)
  if not couple then return "" end
  local kind = couple._virtual and ("v:" .. tostring(couple._viaY or "")) or "r"
  return table.concat({
    tostring(couple.ox or ""),
    tostring(couple.red or ""),
    kind,
    fmt(couple.E0 or 0),
    fmt(couple.n or 0),
    fmt(couple.nH or 0),
    fmt(couple.nH_right or 0)
  }, "/")
end

local function redox_frontier_key(high, low, couple)
  return tostring(high) .. "|" .. tostring(low) .. "|" .. couple_signature(couple)
end

local function has_redox_frontier_at(elementData, pH, cTrace, activeSpeciesSet, convention, key)
  local ladder = M.buildLadder(elementData, pH, cTrace, activeSpeciesSet, convention)
  for i = 2, #ladder do
    local k = redox_frontier_key(ladder[i].form, ladder[i-1].form, ladder[i].couple)
    if k == key then return true end
  end
  return false
end

local function refine_redox_start(elementData, cTrace, activeSpeciesSet, convention, key, firstPH, pHmin, dpH)
  local hi = firstPH
  local lo = math.max(pHmin, firstPH - dpH)
  if lo <= pHmin + 1e-12 then return pHmin end
  if has_redox_frontier_at(elementData, lo, cTrace, activeSpeciesSet, convention, key) then return lo end
  for _ = 1, 40 do
    local mid = 0.5 * (lo + hi)
    if has_redox_frontier_at(elementData, mid, cTrace, activeSpeciesSet, convention, key) then
      hi = mid
    else
      lo = mid
    end
  end
  return hi
end

local function refine_redox_end(elementData, cTrace, activeSpeciesSet, convention, key, lastPH, pHmax, dpH)
  local lo = lastPH
  local hi = math.min(pHmax, lastPH + dpH)
  if hi >= pHmax - 1e-12 then return pHmax end
  if has_redox_frontier_at(elementData, hi, cTrace, activeSpeciesSet, convention, key) then return hi end
  for _ = 1, 40 do
    local mid = 0.5 * (lo + hi)
    if has_redox_frontier_at(elementData, mid, cTrace, activeSpeciesSet, convention, key) then
      lo = mid
    else
      hi = mid
    end
  end
  return lo
end


local function redox_threshold_at(elementData, couple, pH, cTrace, convention, fallback)
  if not couple then return fallback end
  local ox = elementData.species[couple.ox]
  local red = elementData.species[couple.red]
  if not ox or not red then return fallback end
  local E = M.Eeq(couple, ox, red, pH, cTrace, convention)
  if finite(E) then return E end
  return fallback
end

local function ladder_intervals_at(elementData, pHForForms, pHForValues, cTrace, activeSpeciesSet, convention)
  local ladder = M.buildLadder(elementData, pHForForms, cTrace, activeSpeciesSet, convention)
  local intervals = {}
  for i, item in ipairs(ladder) do
    local low = -math.huge
    local high = math.huge
    if i > 1 then
      low = redox_threshold_at(elementData, item.couple, pHForValues, cTrace, convention, item.E_to_prev)
    end
    if i < #ladder then
      local nxt = ladder[i+1]
      high = redox_threshold_at(elementData, nxt.couple, pHForValues, cTrace, convention, nxt.E_to_prev)
    end
    if low and high then
      -- The ladder is ordered at pHForForms.  When pHForValues is extremely
      -- close but not identical, roundoff may swap two almost equal values.
      if high < low and math.abs(high - low) < 1e-9 then
        local mid = 0.5 * (low + high)
        low, high = mid, mid
      end
      intervals[#intervals+1] = {
        form = item.form,
        eMin = low,
        eMax = high,
        ladderIndex = i,
      }
    end
  end
  return intervals
end

local function vertical_frontier_segments(elementData, eq, pHb, cTrace, activeSpeciesSet, range, dpH, convention)
  -- A vertical acid/base boundary is visible only where the dominant species
  -- immediately to the left of the boundary is the acid form and the dominant
  -- species immediately to the right is the base form.  Earlier versions used
  -- the redox ladder exactly at pH = pHb; at a degeneracy, this can assign the
  -- whole vertical extent of an oxidation state and thus draw the boundary
  -- below/above the redox intersection.  We therefore compute the limiting
  -- domains on the two sides of pHb and intersect their E-intervals, while
  -- evaluating the bounding redox potentials at pHb itself.
  local eps = math.min(math.max((dpH or 0.02) * 1e-4, 1e-7), 1e-4)
  local left = ladder_intervals_at(elementData, pHb - eps, pHb, cTrace, activeSpeciesSet, convention)
  local right = ladder_intervals_at(elementData, pHb + eps, pHb, cTrace, activeSpeciesSet, convention)
  local segments = {}
  local function add_intersections(leftId, rightId)
    for _, li in ipairs(left) do
      if li.form == leftId then
        for _, ri in ipairs(right) do
          if ri.form == rightId then
            local e1 = math.max(range.Emin, li.eMin, ri.eMin)
            local e2 = math.min(range.Emax, li.eMax, ri.eMax)
            if finite(e1) and finite(e2) and e2 > e1 + 1e-10 then
              segments[#segments+1] = { eMin = e1, eMax = e2 }
            end
          end
        end
      end
    end
  end
  add_intersections(eq.sideAcid, eq.sideBase)
  -- Robustness for user data with a reversed convention or a merged chain.
  if #segments == 0 then add_intersections(eq.sideBase, eq.sideAcid) end
  table.sort(segments, function(a,b) return a.eMin < b.eMin end)
  return segments
end

local function clip_segment_to_rect(x1, y1, x2, y2, xmin, xmax, ymin, ymax)
  if not (finite(x1) and finite(y1) and finite(x2) and finite(y2)) then return nil end
  local dx, dy = x2 - x1, y2 - y1
  local t0, t1 = 0, 1
  local eps = 1e-14
  local function clip(p, q)
    if math.abs(p) < eps then
      return q >= 0
    end
    local r = q / p
    if p < 0 then
      if r > t1 then return false end
      if r > t0 then t0 = r end
    else
      if r < t0 then return false end
      if r < t1 then t1 = r end
    end
    return true
  end

  if not clip(-dx, x1 - xmin) then return nil end
  if not clip( dx, xmax - x1) then return nil end
  if not clip(-dy, y1 - ymin) then return nil end
  if not clip( dy, ymax - y1) then return nil end
  if t1 < t0 then return nil end

  local ax, ay = x1 + t0 * dx, y1 + t0 * dy
  local bx, by = x1 + t1 * dx, y1 + t1 * dy
  if math.abs(ax - bx) < 1e-12 and math.abs(ay - by) < 1e-12 then return nil end
  return {{ax, ay}, {bx, by}}
end

function M.computeFrontiers(elementData, cTrace, activeSpeciesSet, range, dpH, convention)
  M.resolve_auto_redox_E0s(elementData, elementData and elementData.metal, cTrace, convention)
  convention = M.normalizeConcentrationConvention(convention)
  dpH = dpH or 0.02
  range = M.default_range(elementData, range)
  local pHmin, pHmax = range.pHmin, range.pHmax
  local eMin, eMax = range.Emin, range.Emax

  local present = {}
  for id, _ in pairs(activeSpeciesSet) do
    local sp = elementData.species[id]
    if sp then present[sp.oxState] = true end
  end
  local oxStates = {}
  for n, _ in pairs(present) do oxStates[#oxStates+1] = n end
  table.sort(oxStates)

  -- First pass: detect the pH intervals where each redox boundary exists
  -- chemically. Do not filter by the E window here: chemical existence and
  -- geometric visibility are separate questions. This avoids truncating a
  -- line before its true intersection with the plotting box.
  local segments = {}
  local steps = math.floor((pHmax - pHmin) / dpH + 0.5)
  for k = 0, steps do
    local pH = pHmin + k * dpH
    if pH > pHmax then pH = pHmax end
    local ladder = M.buildLadder(elementData, pH, cTrace, activeSpeciesSet, convention)
    for i = 2, #ladder do
      local high = ladder[i].form
      local low = ladder[i-1].form
      local couple = ladder[i].couple
      local E = ladder[i].E_to_prev
      if finite(E) then
        local key = redox_frontier_key(high, low, couple)
        if not segments[key] then
          segments[key] = { key = key, ox = high, red = low, couple = couple, points = {} }
        end
        local pts = segments[key].points
        pts[#pts+1] = {pH, E}
      end
    end
  end

  local redoxFrontiers = {}
  for _, s in pairs(segments) do
    table.sort(s.points, function(a,b) return a[1] < b[1] end)
    for _, g in ipairs(split_point_groups(s.points, dpH)) do
      local oxSpecies = elementData.species[s.ox]
      local redSpecies = elementData.species[s.red]
      local slope, intercept = nil, nil
      if s.couple and oxSpecies and redSpecies then
        slope, intercept = M.redoxLineCoefficients(s.couple, oxSpecies, redSpecies, cTrace, convention)
      elseif #g >= 2 then
        local p1, p2 = g[1], g[#g]
        if math.abs(p2[1] - p1[1]) > 1e-12 then
          slope = (p2[2] - p1[2]) / (p2[1] - p1[1])
          intercept = p1[2] - slope * p1[1]
        end
      end

      if slope and intercept then
        local pHa = refine_redox_start(elementData, cTrace, activeSpeciesSet, convention, s.key, g[1][1], pHmin, dpH)
        local pHb = refine_redox_end(elementData, cTrace, activeSpeciesSet, convention, s.key, g[#g][1], pHmax, dpH)
        if pHb > pHa + 1e-10 then
          local E1 = slope * pHa + intercept
          local E2 = slope * pHb + intercept
          local visible = clip_segment_to_rect(pHa, E1, pHb, E2, pHmin, pHmax, eMin, eMax)
          if visible then
            redoxFrontiers[#redoxFrontiers+1] = {
              ox = s.ox, red = s.red, couple = s.couple, points = visible,
              fullPoints = {{pHa, E1}, {pHb, E2}},
              slope = slope, intercept = intercept
            }
          end
        end
      end
    end
  end
  table.sort(redoxFrontiers, function(a,b)
    local ax = a.points[1][1]; local bx = b.points[1][1]
    if ax == bx then return (a.ox .. a.red) < (b.ox .. b.red) end
    return ax < bx
  end)

  local acidBaseFrontiers = {}
  for _, n in ipairs(oxStates) do
    local eqs = M.reduceAcidBaseChain(elementData, n, cTrace, activeSpeciesSet, convention)
    for _, eq in ipairs(eqs) do
      local pHb = eq.pHb
      if pHb >= pHmin and pHb <= pHmax then
        local vsegs = vertical_frontier_segments(elementData, eq, pHb, cTrace, activeSpeciesSet, range, dpH, convention)
        for _, vs in ipairs(vsegs) do
          acidBaseFrontiers[#acidBaseFrontiers+1] = {
            acid = eq.sideAcid,
            base = eq.sideBase,
            pH = pHb,
            eMin = vs.eMin,
            eMax = vs.eMax,
            eq = eq,
          }
        end
      end
    end
  end

  table.sort(acidBaseFrontiers, function(a,b)
    if a.pH == b.pH then return (a.acid .. a.base) < (b.acid .. b.base) end
    return a.pH < b.pH
  end)

  return { redoxFrontiers = redoxFrontiers, acidBaseFrontiers = acidBaseFrontiers, range = range }
end

function M.waterFrontiers(pO2, pH2)
  return {
    O2 = function(pH) return 1.23 - 0.06 * pH + 0.06 / 4 * log10(pO2) end,
    H2 = function(pH) return -0.06 * pH - 0.06 / 2 * log10(pH2) end,
  }
end

function M.mhchem_formula(formula)
  local s = formula or "?"
  s = s:gsub("_(%d+)", "%1")
  s = s:gsub("_%{(%d+)%}", "%1")
  return s
end

local function normalize_phase_policy(value)
  local v = trim(tostring(value or "neutral")):lower()
  v = v:gsub("[_-]", " "):gsub("%s+", " ")
  if v == "" or v == "default" or v == "neutral" or v == "neutral only" or v == "only neutral" then
    return "neutral"
  end
  if v == "true" or v == "yes" or v == "on" or v == "1" or v == "all" or v == "always" then
    return "all"
  end
  if v == "false" or v == "no" or v == "off" or v == "0" or v == "none" or v == "never" then
    return "none"
  end
  return "neutral"
end

local function species_is_neutral_label(elementData, id)
  local sp = elementData.species[id]
  if not sp then return true end
  local f = tostring(sp.formula or id or "")
  -- Built-in species formulas encode ionic charge with + or -.
  -- Neutral species such as Fe, Fe(OH)_2, H_2S, V_2O_5 contain neither.
  return not f:find("[+%-]")
end

function M.label_text(elementData, id, phasePolicy)
  local sp = elementData.species[id]
  if not sp then return "?" end
  local f = M.mhchem_formula(sp.formula or id)
  local policy = normalize_phase_policy(phasePolicy)
  local showPhase = false
  if policy == "all" then
    showPhase = true
  elseif policy == "neutral" then
    showPhase = species_is_neutral_label(elementData, id)
  end
  if showPhase and sp.phase and sp.phase ~= "" then
    return "\\ce{" .. f .. "}\\,{\\scriptsize(" .. sp.phase .. ")}"
  end
  return "\\ce{" .. f .. "}"
end
local function canonical_domain_id_text(s)
  s = trim(tostring(s or ""))
  s = s:gsub("^\\ce%s*{%s*(.-)%s*}$", "%1")
  s = s:gsub("%s*%((aq)%)%s*$", "")
  s = s:gsub("%s*%((s)%)%s*$", "")
  s = s:gsub("%s*%((l)%)%s*$", "")
  s = s:gsub("%s*%((g)%)%s*$", "")
  s = s:gsub("_", "")
  s = s:gsub("%^{([^}]*)}", "%1")
  s = s:gsub("%^(%d*[+%-])", "%1")
  s = s:gsub("[{}]", "")
  s = s:gsub("%s+", "")
  return s
end

local function split_selected_domain_list(s)
  local out, cur = {}, {}
  local paren, brace = 0, 0
  s = tostring(s or "")
  for i = 1, #s do
    local ch = s:sub(i,i)
    if ch == "(" then paren = paren + 1
    elseif ch == ")" and paren > 0 then paren = paren - 1
    elseif ch == "{" then brace = brace + 1
    elseif ch == "}" and brace > 0 then brace = brace - 1
    elseif ch == "," and paren == 0 and brace == 0 then
      out[#out+1] = table.concat(cur)
      cur = {}
      goto continue
    end
    cur[#cur+1] = ch
    ::continue::
  end
  local last = table.concat(cur)
  if trim(last) ~= "" then out[#out+1] = last end
  return out
end

local function selected_domain_set(elementData, elementSymbol, list)
  local set = {}
  local raw = trim(list or "")
  if raw == "" then return set end

  local lookup = {}
  for id, sp in pairs(elementData.species or {}) do
    lookup[trim(id)] = id
    lookup[canonical_domain_id_text(id)] = id
    lookup[canonical_domain_id_text(sp.formula or id)] = id
    if sp.phase and sp.phase ~= "" then
      lookup[canonical_domain_id_text((sp.formula or id) .. "(" .. sp.phase .. ")")] = id
      lookup[canonical_domain_id_text(id .. "(" .. sp.phase .. ")")] = id
    end
  end

  for _, item in ipairs(split_selected_domain_list(raw)) do
    local x = trim(item)
    if x ~= "" then
      local id = lookup[x] or lookup[canonical_domain_id_text(x)]
      if id then
        set[id] = true
      else
        lua_warn("Domaine selectionne non reconnu : `" .. x .. "'.")
      end
    end
  end
  return set
end

function M.label_positions(elementData, cTrace, activeSet, range, nx, ny, minCells, convention)
  convention = M.normalizeConcentrationConvention(convention)
  nx = nx or 64
  ny = ny or 64
  minCells = minCells or math.max(5, math.floor(nx * ny * 0.004))
  range = M.default_range(elementData, range)
  local pHmin, pHmax, Emin, Emax = range.pHmin, range.pHmax, range.Emin, range.Emax
  local grid = {}
  for i = 1, nx do
    grid[i] = {}
    local pH = pHmin + (i - 0.5) * (pHmax - pHmin) / nx
    for j = 1, ny do
      local E = Emin + (j - 0.5) * (Emax - Emin) / ny
      grid[i][j] = M.dominantSpeciesAt(elementData, pH, E, cTrace, activeSet, convention)
    end
  end
  local seen = {}
  local labels = {}
  local dirs = {{1,0},{-1,0},{0,1},{0,-1}}
  for i = 1, nx do
    seen[i] = seen[i] or {}
    for j = 1, ny do
      if not seen[i][j] and grid[i][j] then
        local id = grid[i][j]
        local q = {{i,j}}
        seen[i][j] = true
        local head = 1
        local count, sx, sy = 0, 0, 0
        while head <= #q do
          local a, b = q[head][1], q[head][2]
          head = head + 1
          count = count + 1
          sx = sx + (pHmin + (a - 0.5) * (pHmax - pHmin) / nx)
          sy = sy + (Emin + (b - 0.5) * (Emax - Emin) / ny)
          for _, d in ipairs(dirs) do
            local aa, bb = a + d[1], b + d[2]
            if aa >= 1 and aa <= nx and bb >= 1 and bb <= ny then
              seen[aa] = seen[aa] or {}
              if not seen[aa][bb] and grid[aa][bb] == id then
                seen[aa][bb] = true
                q[#q+1] = {aa, bb}
              end
            end
          end
        end
        if count >= minCells then
          labels[#labels+1] = { id = id, pH = sx / count, E = sy / count, count = count }
        end
      end
    end
  end
  table.sort(labels, function(a,b)
    if a.pH == b.pH then return a.E < b.E end
    return a.pH < b.pH
  end)
  return labels
end

M._tex_registry = { last = nil, diagrams = {}, entries = {}, byElement = {} }

function M.annotateFrontiers(elementData, cTrace, fr, convention)
  local byName = {}
  for i, s in ipairs(fr.redoxFrontiers or {}) do
    local name = tostring(i)
    local p1 = s.points and s.points[1]
    local p2 = s.points and s.points[#s.points]
    local oxSpecies = elementData.species[s.ox]
    local redSpecies = elementData.species[s.red]
    local slope, intercept = s.slope, s.intercept
    if (not slope or not intercept) and s.couple and oxSpecies and redSpecies then
      slope, intercept = M.redoxLineCoefficients(s.couple, oxSpecies, redSpecies, cTrace, convention)
    elseif (not slope or not intercept) and p1 and p2 and math.abs(p2[1] - p1[1]) > 1e-12 then
      slope = (p2[2] - p1[2]) / (p2[1] - p1[1])
      intercept = p1[2] - slope * p1[1]
    end
    s.name = name
    s.kind = "redox"
    s.slope = slope
    s.intercept = intercept
    s.horizontal = slope and math.abs(slope) < 1e-10 or false
    byName[name] = {
      kind = "redox",
      name = name,
      ox = s.ox,
      red = s.red,
      slope = slope,
      intercept = intercept,
      horizontal = s.horizontal,
      pH1 = p1 and p1[1] or nil,
      E1 = p1 and p1[2] or nil,
      pH2 = p2 and p2[1] or nil,
      E2 = p2 and p2[2] or nil,
    }
  end

  for i, ab in ipairs(fr.acidBaseFrontiers or {}) do
    local name = letter_name(i)
    ab.name = name
    ab.kind = "vertical"
    byName[name] = {
      kind = "vertical",
      name = name,
      acid = ab.acid,
      base = ab.base,
      pH = ab.pH,
      eMin = ab.eMin,
      eMax = ab.eMax,
    }
  end

  return byName
end

local function canonical_registry_id(element, cTrace, convention)
  return tostring(element or "") .. "@c=" .. fmt_sig(cTrace or 0, 12) .. "@" .. M.normalizeConcentrationConvention(convention or "species")
end

local function parse_registry_selector(sel)
  sel = trim(sel or "")
  if sel == "" then return {raw=""} end
  if not sel:find("=") then return {raw=sel, id=sel} end
  local out = {raw=sel}
  for item in sel:gmatch("([^,]+)") do
    local k, v = item:match("^%s*([^=]+)%s*=%s*(.-)%s*$")
    if k then
      k = trim(k):lower():gsub("[_-]", " "):gsub("%s+", " ")
      v = trim(v)
      if k == "element" or k == "el" or k == "symbol" then out.element = v
      elseif k == "c" or k == "concentration" then out.c = number(v, nil)
      elseif k == "id" or k == "name" then out.id = v
      elseif k == "c convention" or k == "concentration convention" or k == "convention" then out.cConvention = M.normalizeConcentrationConvention(v)
      end
    end
  end
  return out
end

function M.storeFrontierRegistry(id, element, cTrace, fr, byName, convention)
  id = trim(id or "")
  local canonical = canonical_registry_id(element, cTrace, convention)
  local entry = {
    id = id ~= "" and id or canonical,
    canonicalId = canonical,
    element = element,
    cTrace = cTrace,
    cConvention = M.normalizeConcentrationConvention(convention or "species"),
    redoxFrontiers = fr.redoxFrontiers or {},
    acidBaseFrontiers = fr.acidBaseFrontiers or {},
    byName = byName or {},
  }
  M._tex_registry.diagrams[entry.id] = entry
  M._tex_registry.diagrams[canonical] = entry
  M._tex_registry.diagrams[element] = entry
  M._tex_registry.byElement[element] = entry
  M._tex_registry.entries[#M._tex_registry.entries+1] = entry
  M._tex_registry.last = entry.id
  return entry
end

local function registry_entry(id)
  local sel = parse_registry_selector(id)
  if sel.raw == "" then
    local last = M._tex_registry.last
    return last and M._tex_registry.diagrams[last] or nil
  end
  if sel.id and M._tex_registry.diagrams[sel.id] then
    return M._tex_registry.diagrams[sel.id]
  end
  if sel.element and sel.c then
    local conv = sel.cConvention
    for i = #M._tex_registry.entries, 1, -1 do
      local e = M._tex_registry.entries[i]
      if e.element == sel.element and math.abs((e.cTrace or 0) - sel.c) <= math.max(1e-14, 1e-10*math.abs(sel.c))
         and (not conv or e.cConvention == conv) then
        return e
      end
    end
  end
  if sel.element and M._tex_registry.byElement[sel.element] then
    return M._tex_registry.byElement[sel.element]
  end
  if sel.raw and M._tex_registry.diagrams[sel.raw] then
    return M._tex_registry.diagrams[sel.raw]
  end
  return nil
end

local function registry_value(entry, name, field)
  if not entry then return nil end
  local f = entry.byName[tostring(name or "")]
  if not f then return nil end
  return f[field]
end

function M.tex_get_value(id, name, field, sig)
  local entry = registry_entry(id)
  local v = registry_value(entry, name, field)
  if v == nil then return "" end
  if type(v) == "number" then return fmt_sig(v, sig or 10) end
  return tostring(v)
end

function M.tex_get_equation(id, name, sig)
  local entry = registry_entry(id)
  if not entry then return "" end
  local f = entry.byName[tostring(name or "")]
  if not f then return "" end
  if f.kind == "vertical" then
    return "\\mathrm{pH}=" .. fmt_sig(f.pH, sig or 6)
  end
  if f.kind == "redox" and f.slope and f.intercept then
    local m = f.slope
    local b = f.intercept
    local s = "E="
    if math.abs(m) >= 1e-12 then
      s = s .. fmt_sig(m, sig or 6) .. "\\,\\mathrm{pH}"
      if b >= 0 then s = s .. "+" .. fmt_sig(b, sig or 6) else s = s .. fmt_sig(b, sig or 6) end
    else
      s = s .. fmt_sig(b, sig or 6)
    end
    return s
  end
  return ""
end

local function append_unique(list, labels, value, label, tol)
  tol = tol or 1e-8
  for _, v in ipairs(list) do
    if math.abs(v - value) <= tol then return end
  end
  list[#list+1] = value
  labels[#labels+1] = label
end

M._axis_tick_cache = M._axis_tick_cache or {
  x_values = {}, x_labels = {},
  y_values = {}, y_labels = {},
}

function M.tex_reset_axis_ticks()
  M._axis_tick_cache = {
    x_values = {}, x_labels = {},
    y_values = {}, y_labels = {},
  }
end

local function merge_ticks(dst_values, dst_labels, src_values, src_labels)
  for i, v in ipairs(src_values) do
    append_unique(dst_values, dst_labels, v, src_labels[i], 1e-8)
  end
end

local function extra_tick_style_option(axis, style)
  style = trim(style or "")
  if style == "" then return "" end
  -- The option is forwarded to pgfplots' extra tick style for the current
  -- axis.  Typical values are label-oriented TikZ options such as
  -- rotate=90, anchor=east, font=\scriptsize, text=red.
  return string.format(", extra %s tick style={tick label style={%s}}", axis, style)
end

local function emit_extra_ticks(axis, values, labels, style)
  local vtex, ltex = {}, {}
  for i, v in ipairs(values) do
    vtex[#vtex+1] = fmt(v)
    ltex[#ltex+1] = labels[i]
  end
  if axis == "x" then
    tex.print("\\pgfplotsset{extra x ticks={" .. csv(vtex) .. "}, extra x tick labels={" .. csv(ltex) .. "}" .. extra_tick_style_option("x", style) .. "}")
  elseif axis == "y" then
    tex.print("\\pgfplotsset{extra y ticks={" .. csv(vtex) .. "}, extra y tick labels={" .. csv(ltex) .. "}" .. extra_tick_style_option("y", style) .. "}")
  end
end

local function emit_current_axis_extra_ticks(new_x_values, new_x_labels, new_y_values, new_y_labels, x_style, y_style)
  M._axis_tick_cache = M._axis_tick_cache or {x_values={}, x_labels={}, y_values={}, y_labels={}}
  merge_ticks(M._axis_tick_cache.x_values, M._axis_tick_cache.x_labels, new_x_values or {}, new_x_labels or {})
  merge_ticks(M._axis_tick_cache.y_values, M._axis_tick_cache.y_labels, new_y_values or {}, new_y_labels or {})
  emit_extra_ticks("x", M._axis_tick_cache.x_values, M._axis_tick_cache.x_labels, x_style)
  emit_extra_ticks("y", M._axis_tick_cache.y_values, M._axis_tick_cache.y_labels, y_style)
end



-- ---------------------------------------------------------------------------
-- Database extension interface from LaTeX
-- ---------------------------------------------------------------------------

local function tex_escape_message(s)
  s = tostring(s or "")
  s = s:gsub("\\", "\\textbackslash{}")
       :gsub("%%", "\\%%")
       :gsub("#", "\\#")
       :gsub("&", "\\&")
       :gsub("_", "\\_")
       :gsub("%$", "\\$")
       :gsub("{", "\\{")
       :gsub("}", "\\}")
  return s
end

local function lua_warn(msg)
  texio.write_nl("term and log", "pourbaix warning: " .. tostring(msg))
end

local function lua_info(msg)
  texio.write_nl("term and log", "pourbaix: " .. tostring(msg))
end

function M.available_elements()
  local names = {}
  for k, _ in pairs(M.data or {}) do names[#names+1] = k end
  table.sort(names)
  return names
end

function M.available_elements_string()
  return table.concat(M.available_elements(), ", ")
end

local function ensure_element(symbol)
  symbol = trim(symbol or "")
  if symbol == "" then return nil end
  if not M.data[symbol] then
    M.data[symbol] = {
      name = symbol,
      metal = symbol,
      species = {},
      acidBaseEquilibria = {},
      redoxCouples = {},
    }
  else
    M.data[symbol].species = M.data[symbol].species or {}
    M.data[symbol].acidBaseEquilibria = M.data[symbol].acidBaseEquilibria or {}
    M.data[symbol].redoxCouples = M.data[symbol].redoxCouples or {}
  end
  return M.data[symbol]
end

local function parse_optional_number(v)
  if v == nil then return nil end
  v = trim(v)
  if v == "" or v == "auto" then return nil end
  return tonumber(v)
end

local function parse_required_number(v, default)
  local n = parse_optional_number(v)
  if n == nil then return default end
  return n
end

function M.declare_element(spec)
  spec = spec or {}
  local symbol = trim(spec.symbol or spec.element or spec.id or "")
  if symbol == "" then
    lua_warn("\\PourbaixDeclareElement ignored a declaration without a `symbol` key.")
    return false
  end
  local element = ensure_element(symbol)
  element.name = trim(spec.name or element.name or symbol)
  element.metal = trim(spec.metal or element.metal or symbol)
  local Emin = parse_optional_number(spec.Emin or spec.eMin or spec["E min"])
  local Emax = parse_optional_number(spec.Emax or spec.eMax or spec["E max"])
  if Emin and Emax then element.eRange = {Emin, Emax} end
  return true
end

function M.add_species(elementSymbol, spec)
  spec = spec or {}
  elementSymbol = trim(elementSymbol or spec.element or "")
  local element = ensure_element(elementSymbol)
  if not element then
    lua_warn("\\PourbaixAddSpecies ignored a species without an `element` key.")
    return false
  end
  local id = trim(spec.id or "")
  if id == "" then
    lua_warn("\\PourbaixAddSpecies ignored a species without an `id` key for element " .. elementSymbol .. ".")
    return false
  end
  local phase = trim(spec.phase or "aq")
  local oxState = parse_optional_number(spec.oxState or spec.oxstate or spec["ox state"])
  if oxState == nil then
    lua_warn("Species " .. id .. " : missing or non-numeric `ox state`; value 0 used.")
    oxState = 0
  end
  local sp = {
    formula = trim(spec.formula or id),
    phase = phase,
    oxState = oxState,
  }
  local color = trim(spec.color or "")
  if color ~= "" then sp.color = color end
  local nMetal = parse_optional_number(spec.nMetal or spec.nmetal or spec["n metal"])
  if nMetal and math.abs(nMetal - 1) > 1e-12 then sp.nMetal = nMetal end
  element.species[id] = sp
  return true
end

function M.add_acid_base(elementSymbol, eq)
  eq = eq or {}
  elementSymbol = trim(elementSymbol or eq.element or "")
  local element = ensure_element(elementSymbol)
  if not element then
    lua_warn("\\PourbaixAddAcidBase ignored an equilibrium without an `element` key.")
    return false
  end
  local acid = trim(eq.acid or eq.sideAcid or "")
  local base = trim(eq.base or eq.sideBase or "")
  if acid == "" or base == "" then
    lua_warn("\\PourbaixAddAcidBase for " .. elementSymbol .. " : keys `acid` and `base` are required.")
    return false
  end
  local raw_pKa = eq.pKa or eq.pka
  local pKa = parse_optional_number(raw_pKa)
  if pKa == nil then
    if raw_pKa ~= nil and trim(raw_pKa) ~= "" and trim(raw_pKa) ~= "auto" then
      lua_warn("\\PourbaixAddAcidBase for " .. elementSymbol .. " : non-numeric pKa (`" .. tostring(raw_pKa) .. "`); value 0 used.")
    end
    pKa = 0
  end
  local item = {
    sideAcid = acid,
    sideBase = base,
    nH = parse_required_number(eq.nH, 1),
    nH2O = parse_required_number(eq.nH2O or eq.nWater or eq["n water"], 0),
    pKa = pKa,
  }
  local nAcid = parse_optional_number(eq.nAcid or eq.nacid or eq["n acid"])
  local nBase = parse_optional_number(eq.nBase or eq.nbase or eq["n base"])
  if nAcid and math.abs(nAcid - 1) > 1e-12 then item.nAcid = nAcid end
  if nBase and math.abs(nBase - 1) > 1e-12 then item.nBase = nBase end
  table.insert(element.acidBaseEquilibria, item)
  return true
end

function M.add_redox(elementSymbol, couple)
  couple = couple or {}
  elementSymbol = trim(elementSymbol or couple.element or "")
  local element = ensure_element(elementSymbol)
  if not element then
    lua_warn("\\PourbaixAddRedox ignored a couple without an `element` key.")
    return false
  end
  local ox = trim(couple.ox or "")
  local red = trim(couple.red or "")
  if ox == "" or red == "" then
    lua_warn("\\PourbaixAddRedox for " .. elementSymbol .. " : keys `ox` and `red` are required.")
    return false
  end
  local raw_E0 = couple.E0 or couple.e0 or couple.Estandard or couple["E standard"]
  local raw_E0_trim = trim(raw_E0 or "")
  local auto_E0 = (raw_E0_trim == "auto" or raw_E0_trim == "TODO_DERIVED")
  local E0 = parse_optional_number(raw_E0)
  if E0 == nil then
    if raw_E0 ~= nil and raw_E0_trim ~= "" and not auto_E0 then
      lua_warn("\\PourbaixAddRedox for " .. elementSymbol .. " : non-numeric E0 (`" .. tostring(raw_E0) .. "`); value 0 used. Use a numeric value, or `E0={auto}` for a derived potential.")
    end
    E0 = 0
  end
  local item = {
    ox = ox,
    red = red,
    E0 = E0,
    _autoE0 = auto_E0,
    _rawE0 = raw_E0,
    n = parse_required_number(couple.n, 1),
    nH = parse_required_number(couple.nH, 0),
    nH2O = parse_required_number(couple.nH2O or couple.nWater or couple["n water"], 0),
  }
  local fields = {
    {"nOx", "nOx", "nox", "n ox"},
    {"nRed", "nRed", "nred", "n red"},
    {"nH_right", "nH_right", "nHright", "nH right", "n h right"},
    {"nH2O_left", "nH2O_left", "nH2Oleft", "nH2O left", "n water left"},
  }
  for _, names in ipairs(fields) do
    local out = names[1]
    local val = nil
    for i = 2, #names do
      val = couple[names[i]]
      if val ~= nil and trim(val) ~= "" then break end
      val = nil
    end
    local n = parse_optional_number(val)
    if n ~= nil then item[out] = n end
  end
  table.insert(element.redoxCouples, item)
  return true
end

local function merge_element(dst, src)
  dst.name = src.name or dst.name
  dst.metal = src.metal or dst.metal
  dst.eRange = src.eRange or dst.eRange
  dst.species = dst.species or {}
  for id, sp in pairs(src.species or {}) do dst.species[id] = sp end
  dst.acidBaseEquilibria = dst.acidBaseEquilibria or {}
  for _, eq in ipairs(src.acidBaseEquilibria or {}) do table.insert(dst.acidBaseEquilibria, eq) end
  dst.redoxCouples = dst.redoxCouples or {}
  for _, c in ipairs(src.redoxCouples or {}) do table.insert(dst.redoxCouples, c) end
end

function M.merge_data(tbl)
  if type(tbl) ~= "table" then return false end
  for symbol, element in pairs(tbl) do
    if type(element) == "table" then
      local dst = ensure_element(symbol)
      merge_element(dst, element)
    end
  end
  return true
end

function M.load_user_data(filename)
  filename = trim(filename or "")
  if filename == "" then return false end
  local found = kpse.find_file(filename, "tex") or kpse.find_file(filename, "lua") or filename
  local chunk, err = loadfile(found)
  if not chunk then
    lua_warn("Unable to load Lua data file `" .. filename .. "` : " .. tostring(err))
    return false
  end
  local ok, result = pcall(chunk)
  if not ok then
    lua_warn("Error in Lua data file `" .. filename .. "` : " .. tostring(result))
    return false
  end
  if type(result) == "function" then
    ok, result = pcall(result, M)
    if not ok then
      lua_warn("Error while executing the function returned by `" .. filename .. "` : " .. tostring(result))
      return false
    end
  end
  if type(result) == "table" then
    return M.merge_data(result)
  end
  -- A file may also modify require('pourbaix').data directly.
  return true
end

local allowed_phases = { aq = true, s = true, l = true, g = true }

function M.validate_element(symbol)
  symbol = trim(symbol or "")
  local messages = {}
  local element = M.data[symbol]
  if not element then
    return false, {"Element inconnu : " .. symbol .. ". Elements disponibles : " .. M.available_elements_string()}
  end
  element.species = element.species or {}
  element.acidBaseEquilibria = element.acidBaseEquilibria or {}
  element.redoxCouples = element.redoxCouples or {}

  local nSpecies = 0
  for id, sp in pairs(element.species) do
    nSpecies = nSpecies + 1
    if trim(sp.formula or "") == "" then messages[#messages+1] = "Species " .. id .. " : formule absente." end
    if not allowed_phases[trim(sp.phase or "")] then messages[#messages+1] = "Species " .. id .. " : phase inconnue `" .. tostring(sp.phase) .. "`." end
    if tonumber(sp.oxState) == nil then messages[#messages+1] = "Species " .. id .. " : oxState is not numeric." end
    if sp.nMetal ~= nil and (tonumber(sp.nMetal) == nil or tonumber(sp.nMetal) <= 0) then messages[#messages+1] = "Species " .. id .. " : nMetal doit etre strictement positif." end
  end
  if nSpecies == 0 then messages[#messages+1] = "No species defined." end

  for i, eq in ipairs(element.acidBaseEquilibria) do
    local prefix = "Equilibre acido-basique #" .. i .. " : "
    local a, b = element.species[eq.sideAcid], element.species[eq.sideBase]
    if not a then messages[#messages+1] = prefix .. "unknown acid species `" .. tostring(eq.sideAcid) .. "`." end
    if not b then messages[#messages+1] = prefix .. "unknown base species `" .. tostring(eq.sideBase) .. "`." end
    if a and b and tonumber(a.oxState) ~= tonumber(b.oxState) then messages[#messages+1] = prefix .. "the two species do not have the same oxidation state." end
    if tonumber(eq.nH) == nil or tonumber(eq.nH) == 0 then messages[#messages+1] = prefix .. "nH doit etre non nul." end
    if tonumber(eq.pKa) == nil then messages[#messages+1] = prefix .. "pKa is not numeric." end
  end

  for i, c in ipairs(element.redoxCouples) do
    local prefix = "Couple redox #" .. i .. " : "
    local ox, red = element.species[c.ox], element.species[c.red]
    if not ox then messages[#messages+1] = prefix .. "unknown oxidized species `" .. tostring(c.ox) .. "`." end
    if not red then messages[#messages+1] = prefix .. "unknown reduced species `" .. tostring(c.red) .. "`." end
    if ox and red and tonumber(ox.oxState) and tonumber(red.oxState) and tonumber(ox.oxState) <= tonumber(red.oxState) then
      messages[#messages+1] = prefix .. "oxState(ox) doit etre strictement superieur a oxState(red)."
    end
    if tonumber(c.E0) == nil then messages[#messages+1] = prefix .. "E0 is not numeric." end
    if tonumber(c.n) == nil or tonumber(c.n) <= 0 then messages[#messages+1] = prefix .. "n doit etre strictement positif." end
  end

  return #messages == 0, messages
end

function M.tex_validate_element(symbol)
  local ok, messages = M.validate_element(symbol)
  if ok then
    lua_info("Validation de l'element " .. tostring(symbol) .. " : OK.")
  else
    for _, msg in ipairs(messages) do lua_warn(msg) end
  end
end

local function dump_value(v, indent, seen)
  indent = indent or ""
  seen = seen or {}
  if type(v) ~= "table" then
    if type(v) == "string" then return string.format("%q", v) end
    return tostring(v)
  end
  if seen[v] then return "<cycle>" end
  seen[v] = true
  local parts = {"{"}
  local keys = {}
  for k, _ in pairs(v) do keys[#keys+1] = k end
  table.sort(keys, function(a,b) return tostring(a) < tostring(b) end)
  for _, k in ipairs(keys) do
    parts[#parts+1] = indent .. "  [" .. dump_value(k, "", seen) .. "] = " .. dump_value(v[k], indent .. "  ", seen) .. ","
  end
  parts[#parts+1] = indent .. "}"
  seen[v] = nil
  return table.concat(parts, "\n")
end

function M.tex_show_element_data(symbol)
  symbol = trim(symbol or "")
  local element = M.data[symbol]
  if not element then
    lua_warn("Element inconnu : " .. symbol)
    return
  end
  texio.write_nl("term and log", "--- Normalized Pourbaix data for " .. symbol .. " ---")
  texio.write_nl("term and log", dump_value(element, ""))
  texio.write_nl("term and log", "--- End of Pourbaix data for " .. symbol .. " ---")
end


-- ---------------------------------------------------------------------------
-- Thermodynamic skeleton generator from a species list.
-- This function does not modify the database: it writes commands to the log
-- for the user to complete: \PourbaixAddSpecies, \PourbaixAddAcidBase,
-- and \PourbaixAddRedox.
-- ---------------------------------------------------------------------------

local function gcd_int(a, b)
  a, b = math.abs(math.floor(a or 0)), math.abs(math.floor(b or 0))
  while b ~= 0 do a, b = b, a % b end
  if a == 0 then return 1 end
  return a
end

local function lcm_int(a, b)
  a, b = math.abs(math.floor(a or 1)), math.abs(math.floor(b or 1))
  return math.floor(a / gcd_int(a, b) * b)
end

local function split_top_level_commas(s)
  s = tostring(s or "")
  local out, buf = {}, {}
  local par, bra = 0, 0
  for i = 1, #s do
    local ch = s:sub(i,i)
    if ch == "(" then par = par + 1
    elseif ch == ")" and par > 0 then par = par - 1
    elseif ch == "{" then bra = bra + 1
    elseif ch == "}" and bra > 0 then bra = bra - 1
    elseif ch == "," and par == 0 and bra == 0 then
      local item = trim(table.concat(buf))
      if item ~= "" then out[#out+1] = item end
      buf = {}
      goto continue
    end
    buf[#buf+1] = ch
    ::continue::
  end
  local item = trim(table.concat(buf))
  if item ~= "" then out[#out+1] = item end
  return out
end

local function parse_charge(raw)
  raw = trim(raw or "")
  if raw == "" then return 0, "" end
  local sign = raw:sub(-1)
  if sign ~= "+" and sign ~= "-" then return 0, "" end
  local mag = raw:sub(1, -2)
  if mag == "" then mag = 1 else mag = tonumber(mag) or 1 end
  local val = (sign == "+") and mag or -mag
  local label = (math.abs(val) == 1) and sign or (tostring(math.abs(val)) .. sign)
  return val, label
end

local function strip_tex_formula(s)
  s = trim(s or "")
  s = s:gsub("^\\ce%s*{%s*(.-)%s*}$", "%1")
  s = s:gsub("%s+", "")
  s = s:gsub("_%{?([0-9]+)%}?", "%1")
  s = s:gsub("%^{([0-9]*[+%-])}", "^%1")
  return s
end

local function split_phase_and_charge(raw, default_phase)
  local s = strip_tex_formula(raw)
  local phase = trim(default_phase or "aq")
  local ph = s:match("%((aq|s|l|g)%)$")
  -- Lua patterns do not support alternation in captures; handle manually.
  ph = s:match("%((aq)%)$") or s:match("%((s)%)$") or s:match("%((l)%)$") or s:match("%((g)%)$")
  if ph then
    phase = ph
    s = s:gsub("%(" .. ph .. "%)$", "")
  end

  local core, charge_raw = s:match("^(.-)%^%{?([0-9]*[+%-])%}?$")
  if not core then
    -- Also accepts a terminal charge without a caret, e.g. Fe2+ or ClO-.
    local c2, q2 = s:match("^(.-)([0-9]*[+%-])$")
    if c2 and c2 ~= "" then core, charge_raw = c2, q2 end
  end
  if not core then core, charge_raw = s, "" end
  local charge, charge_label = parse_charge(charge_raw)
  return core, phase, charge, charge_label
end

local function parse_atom_counts(formula)
  formula = strip_tex_formula(formula or "")
  local stack = {{} }
  local i = 1
  local function add(el, n)
    local top = stack[#stack]
    top[el] = (top[el] or 0) + n
  end
  local function parse_number()
    local j = i
    while j <= #formula and formula:sub(j,j):match("%d") do j = j + 1 end
    if j > i then
      local n = tonumber(formula:sub(i, j-1)) or 1
      i = j
      return n
    end
    return 1
  end
  while i <= #formula do
    local ch = formula:sub(i,i)
    if ch == "(" then
      stack[#stack+1] = {}
      i = i + 1
    elseif ch == ")" then
      i = i + 1
      local mult = parse_number()
      local top = table.remove(stack)
      if not top then top = {} end
      for el, n in pairs(top) do add(el, n * mult) end
    elseif ch:match("%u") then
      local el = ch
      i = i + 1
      if i <= #formula and formula:sub(i,i):match("%l") then
        el = el .. formula:sub(i,i)
        i = i + 1
      end
      local n = parse_number()
      add(el, n)
    else
      i = i + 1
    end
  end
  return stack[1] or {}
end

local function infer_oxidation_state(element, counts, charge)
  local nu = counts[element] or 0
  if nu <= 0 then return nil end
  local known = 0
  for el, n in pairs(counts) do
    if el ~= element then
      if el == "O" then known = known - 2*n
      elseif el == "H" then known = known + n
      else known = known + 0 end
    end
  end
  local x = (charge - known) / nu
  if math.abs(x - math.floor(x + 0.5)) < 1e-10 then x = math.floor(x + 0.5) end
  return x
end

local function latex_formula_from_core(core, charge_label)
  local out = {}
  local i = 1
  while i <= #core do
    local ch = core:sub(i,i)
    if ch:match("%d") then
      local j = i
      while j <= #core and core:sub(j,j):match("%d") do j = j + 1 end
      local n = core:sub(i, j-1)
      out[#out+1] = (n:len() == 1) and ("_" .. n) or ("_{" .. n .. "}")
      i = j
    else
      out[#out+1] = ch
      i = i + 1
    end
  end
  local f = table.concat(out)
  if charge_label and charge_label ~= "" then
    f = f .. "^{" .. charge_label .. "}"
  end
  return f
end

local function id_from_core(element, core, charge_label)
  local id = core
  if charge_label and charge_label ~= "" then
    local simple_element = (core == element)
    local needs_space = (not simple_element) and charge_label:match("^%d")
    if needs_space then id = id .. " " .. charge_label else id = id .. charge_label end
  end
  return id
end

local function parse_template_species(raw, element, default_phase)
  local core, phase, charge, charge_label = split_phase_and_charge(raw, default_phase)
  local counts = parse_atom_counts(core)
  local nMetal = counts[element] or 0
  local ox = infer_oxidation_state(element, counts, charge)
  local id = id_from_core(element, core, charge_label)
  return {
    raw = raw,
    core = core,
    phase = phase,
    charge = charge,
    chargeLabel = charge_label,
    counts = counts,
    nMetal = nMetal > 0 and nMetal or 1,
    oxState = ox,
    id = id,
    formula = latex_formula_from_core(core, charge_label),
  }
end

local function infer_element_from_species(items)
  for _, raw in ipairs(items) do
    local core = split_phase_and_charge(raw, "aq")
    local counts = parse_atom_counts(core)
    local best, bestn = nil, 0
    for el, n in pairs(counts) do
      if el ~= "H" and el ~= "O" and n > bestn then best, bestn = el, n end
    end
    if best then return best end
  end
  return "X"
end

local function fmt_int_or_num(x)
  if x == nil then return "TODO" end
  if math.abs(x - math.floor(x + 0.5)) < 1e-10 then return tostring(math.floor(x + 0.5)) end
  return fmt(x)
end

local function command_quote_id(id)
  if tostring(id):match("^[A-Za-z][A-Za-z0-9]*$") then return tostring(id) end
  return "{" .. tostring(id) .. "}"
end

local function template_line(line)
  texio.write_nl("term and log", line)
end

local function species_label(sp)
  return sp.id or sp.raw or "?"
end

local function species_equation_label(sp, with_phase)
  if not sp then return "?" end
  local core = sp.core or sp.id or sp.raw or "?"
  local charge = sp.chargeLabel or ""
  local out = core
  if charge ~= "" then out = out .. "^" .. charge end
  if with_phase and sp.phase and sp.phase ~= "" then out = out .. "(" .. sp.phase .. ")" end
  return out
end

local function term_for_equation(coef, label)
  coef = coef or 1
  if math.abs(coef) < 1e-12 then return nil end
  if math.abs(coef - 1) < 1e-12 then return label end
  return fmt_int_or_num(coef) .. " " .. label
end

local function join_equation_terms(terms)
  local out = {}
  for _,t in ipairs(terms or {}) do if t and t ~= "" then out[#out+1] = t end end
  if #out == 0 then return "0" end
  return table.concat(out, " + ")
end

local function acid_base_equation_text(ab)
  local left = { term_for_equation(ab.nAcid or 1, species_equation_label(ab.acid, true)) }
  if ab.nH2O and math.abs(ab.nH2O) > 1e-12 then left[#left+1] = term_for_equation(ab.nH2O, "H2O") end
  local right = { term_for_equation(ab.nBase or 1, species_equation_label(ab.base, true)) }
  if ab.nH and math.abs(ab.nH) > 1e-12 then right[#right+1] = term_for_equation(ab.nH, "H+") end
  return join_equation_terms(left) .. " <=> " .. join_equation_terms(right)
end

local function redox_equation_text(r)
  local left = { term_for_equation(r.nOx or 1, species_equation_label(r.ox, true)) }
  if r.nH and math.abs(r.nH) > 1e-12 then left[#left+1] = term_for_equation(r.nH, "H+") end
  if r.nH2O_left and math.abs(r.nH2O_left) > 1e-12 then left[#left+1] = term_for_equation(r.nH2O_left, "H2O") end
  if r.n and math.abs(r.n) > 1e-12 then left[#left+1] = term_for_equation(r.n, "e-") end
  local right = { term_for_equation(r.nRed or 1, species_equation_label(r.red, true)) }
  if r.nH_right and math.abs(r.nH_right) > 1e-12 then right[#right+1] = term_for_equation(r.nH_right, "H+") end
  if r.nH2O and math.abs(r.nH2O) > 1e-12 then right[#right+1] = term_for_equation(r.nH2O, "H2O") end
  return join_equation_terms(left) .. " <=> " .. join_equation_terms(right)
end

local function ab_label(ab)
  return "pKa(" .. species_label(ab.acid) .. "/" .. species_label(ab.base) .. ")"
end
local function redox_label(r)
  return "E0(" .. species_label(r.ox) .. "/" .. species_label(r.red) .. ")"
end


local function infer_acid_base(a, b)
  local lcm = lcm_int(a.nMetal, b.nMetal)
  local ca, cb = lcm / a.nMetal, lcm / b.nMetal
  local function try(acid, base, cAcid, cBase)
    local oA = acid.counts.O or 0
    local oB = base.counts.O or 0
    local hA = acid.counts.H or 0
    local hB = base.counts.H or 0
    local w = cBase * oB - cAcid * oA
    local nH = cAcid * hA + 2*w - cBase * hB
    local lhsCharge = cAcid * acid.charge
    local rhsCharge = cBase * base.charge + nH
    if w >= -1e-10 and nH > 1e-10 and math.abs(lhsCharge - rhsCharge) < 1e-8 then
      return {
        acid = acid, base = base,
        nAcid = cAcid, nBase = cBase,
        nH = math.max(0, nH), nH2O = math.max(0, w)
      }
    end
    return nil
  end
  return try(a, b, ca, cb) or try(b, a, cb, ca)
end

local function infer_redox(ox, red)
  local lcm = lcm_int(ox.nMetal, red.nMetal)
  local cOx, cRed = lcm / ox.nMetal, lcm / red.nMetal
  local oOx, oRed = ox.counts.O or 0, red.counts.O or 0
  local hOx, hRed = ox.counts.H or 0, red.counts.H or 0
  local deltaO = cRed * oRed - cOx * oOx
  local nH2Oleft, nH2Oright = 0, 0
  if deltaO > 0 then nH2Oleft = deltaO else nH2Oright = -deltaO end
  local leftH = cOx*hOx + 2*nH2Oleft
  local rightH = cRed*hRed + 2*nH2Oright
  local nHleft, nHright = 0, 0
  if rightH > leftH then nHleft = rightH - leftH else nHright = leftH - rightH end
  local leftCharge = cOx*ox.charge + nHleft
  local rightCharge = cRed*red.charge + nHright
  local ne = leftCharge - rightCharge
  if ne <= 1e-10 then return nil end
  return {
    ox = ox, red = red,
    n = ne, nOx = cOx, nRed = cRed,
    nH = nHleft, nH_right = nHright,
    nH2O = nH2Oright, nH2O_left = nH2Oleft,
  }
end

local function emit_keyval_command(name, parts)
  template_line("\\" .. name .. "[")
  for i, kv in ipairs(parts) do
    local comma = (i < #parts) and "," or ""
    template_line("  " .. kv .. comma)
  end
  template_line("]")
end

local function emit_species_template(element, sp)
  emit_keyval_command("PourbaixAddSpecies", {
    string.format("element=%s", element),
    string.format("id=%s", command_quote_id(sp.id)),
    string.format("formula={%s}", sp.formula),
    string.format("phase=%s", sp.phase),
    string.format("ox state=%s", fmt_int_or_num(sp.oxState)),
    string.format("n metal=%s", fmt_int_or_num(sp.nMetal)),
  })
end

local function emit_acidbase_template(element, ab, pka_text)
  template_line("% equation associated with " .. ab_label(ab) .. " : " .. acid_base_equation_text(ab))
  local parts = {
    string.format("element=%s", element),
    string.format("acid=%s", command_quote_id(ab.acid.id)),
    string.format("base=%s", command_quote_id(ab.base.id)),
  }
  if math.abs(ab.nAcid - 1) > 1e-12 then parts[#parts+1] = "n acid=" .. fmt_int_or_num(ab.nAcid) end
  if math.abs(ab.nBase - 1) > 1e-12 then parts[#parts+1] = "n base=" .. fmt_int_or_num(ab.nBase) end
  parts[#parts+1] = "nH=" .. fmt_int_or_num(ab.nH)
  parts[#parts+1] = "nH2O=" .. fmt_int_or_num(ab.nH2O)
  parts[#parts+1] = "pKa={" .. (pka_text or "TODO") .. "}"
  emit_keyval_command("PourbaixAddAcidBase", parts)
end

local function emit_redox_template(element, r, e0_text, note)
  if note and note ~= "" then template_line("% " .. note) end
  template_line("% half-equation associated with " .. redox_label(r) .. " : " .. redox_equation_text(r))
  local parts = {
    string.format("element=%s", element),
    string.format("ox=%s", command_quote_id(r.ox.id)),
    string.format("red=%s", command_quote_id(r.red.id)),
    "E0={" .. (e0_text or "TODO") .. "}",
    "n=" .. fmt_int_or_num(r.n),
  }
  if math.abs(r.nOx - 1) > 1e-12 then parts[#parts+1] = "n ox=" .. fmt_int_or_num(r.nOx) end
  if math.abs(r.nRed - 1) > 1e-12 then parts[#parts+1] = "n red=" .. fmt_int_or_num(r.nRed) end
  if r.nH and math.abs(r.nH) > 1e-12 then parts[#parts+1] = "nH=" .. fmt_int_or_num(r.nH) end
  if r.nH_right and math.abs(r.nH_right) > 1e-12 then parts[#parts+1] = "nH right=" .. fmt_int_or_num(r.nH_right) end
  if r.nH2O and math.abs(r.nH2O) > 1e-12 then parts[#parts+1] = "nH2O=" .. fmt_int_or_num(r.nH2O) end
  if r.nH2O_left and math.abs(r.nH2O_left) > 1e-12 then parts[#parts+1] = "nH2O left=" .. fmt_int_or_num(r.nH2O_left) end
  emit_keyval_command("PourbaixAddRedox", parts)
end

-- Linear expressions used by the skeleton generator. An expression represents
-- a linear combination of independent constants plus a numerical constant.
-- Coefficients are in volts; pKa values are converted through
-- S_TEMPLATE = 0.06 V to stay consistent with the plotting engine.
local S_TEMPLATE = 0.06

local function expr_zero() return {c = 0, t = {}} end
local function expr_const(c) return {c = c or 0, t = {}} end
local function expr_var(name) return {c = 0, t = {[name] = 1}} end
local function expr_clone(a)
  local r = {c = a.c or 0, t = {}}
  for k,v in pairs(a.t or {}) do r.t[k] = v end
  return r
end
local function expr_add(a,b)
  local r = expr_clone(a or expr_zero())
  r.c = r.c + (b and b.c or 0)
  for k,v in pairs((b and b.t) or {}) do r.t[k] = (r.t[k] or 0) + v end
  return r
end
local function expr_sub(a,b) return expr_add(a, expr_scale(b, -1)) end
function expr_scale(a, s)
  local r = {c = (a and a.c or 0) * s, t = {}}
  for k,v in pairs((a and a.t) or {}) do r.t[k] = v * s end
  return r
end
local function expr_eval(a, values)
  local x = a.c or 0
  for k,v in pairs(a.t or {}) do
    local val = values and values[k]
    if val == nil then return nil end
    x = x + v * val
  end
  return x
end
local function fmt_coef(x)
  if math.abs(x) < 5e-12 then return nil end
  if math.abs(x - 1) < 5e-12 then return "" end
  if math.abs(x + 1) < 5e-12 then return "-" end
  return string.format("%.10g*", x)
end
local function expr_format(a)
  local parts = {}
  local keys = {}
  for k,v in pairs(a.t or {}) do if math.abs(v) > 5e-12 then keys[#keys+1]=k end end
  table.sort(keys)
  for _,k in ipairs(keys) do
    local coef = a.t[k]
    local s = fmt_coef(coef)
    if s then parts[#parts+1] = s .. k end
  end
  if math.abs(a.c or 0) > 5e-12 or #parts == 0 then parts[#parts+1] = string.format("%.10g", a.c or 0) end
  local out = table.concat(parts, " + ")
  out = out:gsub("%+ %-", "- ")
  out = out:gsub("%-([%w])", "- %1")
  return out
end

local function pair_key(a,b) return (a.id or "?") .. "//" .. (b.id or "?") end

local function acid_base_score(sp)
  return 10*(sp.counts.O or 0) + (sp.counts.H or 0) - 0.01*(sp.charge or 0)
end

local function build_acid_base_forest(species)
  local byOx = {}
  for _, sp in ipairs(species) do
    if sp.oxState then
      byOx[sp.oxState] = byOx[sp.oxState] or {}
      byOx[sp.oxState][#byOx[sp.oxState]+1] = sp
    end
  end
  local forest = {groups = {}, edges = {}, edge_by_key = {}}
  local oxLevels = {}
  for ox,_ in pairs(byOx) do oxLevels[#oxLevels+1]=ox end
  table.sort(oxLevels)
  for _,ox in ipairs(oxLevels) do
    local group = byOx[ox]
    table.sort(group, function(a,b)
      local sa, sb = acid_base_score(a), acid_base_score(b)
      if math.abs(sa-sb) < 1e-12 then return a.id < b.id end
      return sa < sb
    end)
    local ginfo = {ox = ox, species = group, root = group[1], rank = {}}
    for i,sp in ipairs(group) do ginfo.rank[sp.id] = i-1 end
    forest.groups[ox] = ginfo
    -- Keep a spanning tree in the acid -> base direction: each new species is
    -- connected to the nearest previous species for which an equilibrium can
    -- be inferred. For classical hydroxide/oxide diagrams this gives the
    -- expected independent vertical constants.
    for i=2,#group do
      local chosen = nil
      for j=i-1,1,-1 do
        local e = infer_acid_base(group[j], group[i])
        if e then chosen = e; break end
      end
      if chosen then
        chosen.independent = true
        chosen.label = ab_label(chosen)
        forest.edges[#forest.edges+1] = chosen
        forest.edge_by_key[pair_key(chosen.acid, chosen.base)] = chosen
      end
    end
  end
  return forest
end

local function add_ab_equation(eqs, ab)
  local rhs = expr_scale(expr_var(ab.label), S_TEMPLATE)
  eqs[#eqs+1] = {
    type="ab", a=ab.acid, b=ab.base, nA=ab.nAcid or 1, nB=ab.nBase or 1,
    rhs=rhs, label=ab.label
  }
end

local function add_redox_equation(eqs, r, label)
  local rhs = expr_scale(expr_var(label), r.n or 1)
  eqs[#eqs+1] = {
    type="redox", ox=r.ox, red=r.red, nOx=r.nOx or 1, nRed=r.nRed or 1,
    rhs=rhs, label=label
  }
end

local function solve_g_expressions(species, ab_edges, independent_redox)
  local g = {}
  local eqs = {}
  for _,ab in ipairs(ab_edges) do add_ab_equation(eqs, ab) end
  for _,it in ipairs(independent_redox) do add_redox_equation(eqs, it.r, it.label) end
  table.sort(species, function(a,b)
    local oa, ob = a.oxState or 0, b.oxState or 0
    if oa == ob then return acid_base_score(a) < acid_base_score(b) end
    return oa < ob
  end)
  if species[1] then g[species[1].id] = expr_zero() end
  local changed = true
  local guard = 0
  while changed and guard < 200 do
    changed = false; guard = guard + 1
    for _,eq in ipairs(eqs) do
      if eq.type == "ab" then
        local ga, gb = g[eq.a.id], g[eq.b.id]
        if ga and not gb then
          -- nB*gB - nA*gA = rhs
          g[eq.b.id] = expr_scale(expr_add(expr_scale(ga, eq.nA), eq.rhs), 1/eq.nB)
          changed = true
        elseif gb and not ga then
          -- nA*gA = nB*gB - rhs
          g[eq.a.id] = expr_scale(expr_sub(expr_scale(gb, eq.nB), eq.rhs), 1/eq.nA)
          changed = true
        end
      else
        local gox, gred = g[eq.ox.id], g[eq.red.id]
        if gred and not gox then
          -- nOx*gOx - nRed*gRed = rhs
          g[eq.ox.id] = expr_scale(expr_add(expr_scale(gred, eq.nRed), eq.rhs), 1/eq.nOx)
          changed = true
        elseif gox and not gred then
          -- nRed*gRed = nOx*gOx - rhs
          g[eq.red.id] = expr_scale(expr_sub(expr_scale(gox, eq.nOx), eq.rhs), 1/eq.nRed)
          changed = true
        end
      end
    end
  end
  return g
end

local function build_independent_redox(species, forest)
  local levels = {}
  for ox,_ in pairs(forest.groups) do levels[#levels+1]=ox end
  table.sort(levels)
  local out = {}
  for i=2,#levels do
    local high = forest.groups[levels[i]].root
    local low  = forest.groups[levels[i-1]].root
    local r = infer_redox(high, low)
    if r then
      local label = redox_label(r)
      out[#out+1] = {r=r, label=label, key=pair_key(r.ox, r.red)}
    end
  end
  return out
end

local function should_emit_redox_pair(mode, forest, high, low, adjacentOx)
  if mode == "all" then return true end
  if mode == "adjacent" then return adjacentOx end
  -- dominant mode: generate only couples between adjacent oxidation states.
  -- Boundaries between non-adjacent states are not entered as explicit data: if
  -- an intermediate state is unstable, buildLadder dynamically reconstructs
  -- it as a virtual disproportionation couple. This prevents stable
  -- intermediates, e.g. V(III) / V2O3 in vanadium, from disappearing.
  return adjacentOx
end

local function find_database_redox_E0(elementSymbol, ox, red)
  local data = M.data[elementSymbol]
  if not data then return nil end
  for _,c in ipairs(data.redoxCouples or {}) do
    if c.ox == ox and c.red == red then return c.E0 end
  end
  return nil
end

local function find_database_pKa(elementSymbol, acid, base)
  local data = M.data[elementSymbol]
  if not data then return nil end
  for _,ab in ipairs(data.acidBaseEquilibria or {}) do
    if ab.sideAcid == acid and ab.sideBase == base then return ab.pKa end
  end
  return nil
end

local function build_runtime_species_objects(elementSymbol, elementData)
  local species = {}
  local metal = trim((elementData and elementData.metal) or elementSymbol or "")
  if metal == "" then metal = elementSymbol end
  for id, sp in pairs(elementData.species or {}) do
    local raw = id
    if sp.phase and sp.phase ~= "" and sp.phase ~= "aq" then raw = raw .. "(" .. sp.phase .. ")" end
    local tsp = parse_template_species(raw, metal, sp.phase or "aq")
    tsp.id = id
    tsp.formula = sp.formula or tsp.formula
    tsp.phase = sp.phase or tsp.phase
    tsp.oxState = sp.oxState or tsp.oxState
    tsp.nMetal = sp.nMetal or tsp.nMetal or 1
    species[#species+1] = tsp
  end
  table.sort(species, function(a,b)
    local oa, ob = a.oxState or 0, b.oxState or 0
    if oa == ob then
      local sa, sb = acid_base_score(a), acid_base_score(b)
      if math.abs(sa-sb) < 1e-12 then return a.id < b.id end
      return sa < sb
    end
    return oa < ob
  end)
  return species
end

local function solve_effective_species_g(elementSymbol, elementData, cTrace, convention)
  -- Unknowns are not the standard potentials of species, but the effective
  -- intercept potentials H_i = G_i + S log a_i associated with the current
  -- plotting convention.  This is crucial for local concentration conventions:
  -- it forces triple points obtained by an acid/base boundary and a redox
  -- boundary to close in the actual diagram.
  local species = build_runtime_species_objects(elementSymbol, elementData)
  if #species == 0 then return nil, nil end
  local byId = {}
  for _,sp in ipairs(species) do byId[sp.id] = sp end
  local eqs = {}
  for _,ab in ipairs(elementData.acidBaseEquilibria or {}) do
    local acid, base = byId[ab.sideAcid], byId[ab.sideBase]
    if acid and base then
      local nA, nB, nH = ab.nAcid or 1, ab.nBase or 1, ab.nH or 1
      local logB = M.logActivity(elementData.species[ab.sideBase], cTrace, convention, elementData.species[ab.sideAcid])
      local logA = M.logActivity(elementData.species[ab.sideAcid], cTrace, convention, elementData.species[ab.sideBase])
      local pHb = ((ab.pKa or 0) + nB * logB - nA * logA) / nH
      eqs[#eqs+1] = { type='ab', a=acid, b=base, nA=nA, nB=nB, rhs=expr_const(S_TEMPLATE * nH * pHb) }
    end
  end
  for _,c in ipairs(elementData.redoxCouples or {}) do
    if not c._autoE0 then
      local ox, red = byId[c.ox], byId[c.red]
      if ox and red then
        local _, intercept = M.redoxLineCoefficients(c, elementData.species[c.ox], elementData.species[c.red], cTrace, convention)
        eqs[#eqs+1] = { type='redox', ox=ox, red=red, nOx=c.nOx or 1, nRed=c.nRed or 1, rhs=expr_const((c.n or 1) * intercept) }
      end
    end
  end
  local g = {}
  if species[1] then g[species[1].id] = expr_zero() end
  local changed = true
  local guard = 0
  while changed and guard < 500 do
    changed = false; guard = guard + 1
    for _,eq in ipairs(eqs) do
      if eq.type == 'ab' then
        local ga, gb = g[eq.a.id], g[eq.b.id]
        if ga and not gb then
          g[eq.b.id] = expr_scale(expr_add(expr_scale(ga, eq.nA), eq.rhs), 1/eq.nB)
          changed = true
        elseif gb and not ga then
          g[eq.a.id] = expr_scale(expr_sub(expr_scale(gb, eq.nB), eq.rhs), 1/eq.nA)
          changed = true
        end
      else
        local gox, gred = g[eq.ox.id], g[eq.red.id]
        if gred and not gox then
          g[eq.ox.id] = expr_scale(expr_add(expr_scale(gred, eq.nRed), eq.rhs), 1/eq.nOx)
          changed = true
        elseif gox and not gred then
          g[eq.red.id] = expr_scale(expr_sub(expr_scale(gox, eq.nOx), eq.rhs), 1/eq.nRed)
          changed = true
        end
      end
    end
  end
  return g, byId
end

local function add_missing_adjacent_auto_couples(elementData, elementSymbol)
  local species = build_runtime_species_objects(elementSymbol or elementData.metal or "X", elementData)
  if #species == 0 then return end
  local levels, byLevel, seen = {}, {}, {}
  for _,sp in ipairs(species) do
    if sp.oxState then
      byLevel[sp.oxState] = byLevel[sp.oxState] or {}
      byLevel[sp.oxState][#byLevel[sp.oxState]+1] = sp
      if not seen[sp.oxState] then levels[#levels+1] = sp.oxState; seen[sp.oxState] = true end
    end
  end
  table.sort(levels)
  local existing = {}
  for _,c in ipairs(elementData.redoxCouples or {}) do existing[tostring(c.ox).."//"..tostring(c.red)] = true end
  for i=2,#levels do
    local highList = byLevel[levels[i]] or {}
    local lowList  = byLevel[levels[i-1]] or {}
    for _,ox in ipairs(highList) do
      for _,red in ipairs(lowList) do
        local key = ox.id .. "//" .. red.id
        if not existing[key] then
          local r = infer_redox(ox, red)
          if r then
            local item = {
              ox = r.ox.id, red = r.red.id, E0 = 0,
              n = r.n, nOx = r.nOx, nRed = r.nRed,
              nH = r.nH, nH_right = r.nH_right,
              nH2O = r.nH2O, nH2O_left = r.nH2O_left,
              _autoE0 = true, _runtimeAuto = true,
            }
            elementData.redoxCouples[#elementData.redoxCouples+1] = item
            existing[key] = true
          end
        end
      end
    end
  end
end

function M.resolve_auto_redox_E0s(elementData, elementSymbol, cTrace, convention)
  if not elementData then return end
  cTrace = cTrace or 1e-2
  convention = M.normalizeConcentrationConvention(convention or "species")
  local hasAuto = false
  for _,c in ipairs(elementData.redoxCouples or {}) do if c._autoE0 then hasAuto = true; break end end
  if not hasAuto then return end
  add_missing_adjacent_auto_couples(elementData, elementSymbol or elementData.metal or "X")
  local g, byId = solve_effective_species_g(elementSymbol or elementData.metal or "X", elementData, cTrace, convention)
  if not g then return end
  for _,c in ipairs(elementData.redoxCouples or {}) do
    if c._autoE0 then
      local ox, red = byId[c.ox], byId[c.red]
      local gox, gred = ox and g[ox.id], red and g[red.id]
      if ox and red and gox and gred then
        local nOx, nRed, n = c.nOx or 1, c.nRed or 1, c.n or 1
        local interceptExpr = expr_scale(expr_sub(expr_scale(gox, nOx), expr_scale(gred, nRed)), 1/n)
        local intercept = expr_eval(interceptExpr, {})
        if intercept ~= nil then
          local logOx = M.logActivity(elementData.species[c.ox], cTrace, convention, elementData.species[c.red])
          local logRed = M.logActivity(elementData.species[c.red], cTrace, convention, elementData.species[c.ox])
          local activityTerm = (S_TEMPLATE / n) * (nOx * logOx - nRed * logRed)
          c.E0 = intercept - activityTerm
          c._autoE0_resolved = true
        else
          lua_warn("Unable to resolve the derived potential automatically for couple " .. tostring(c.ox) .. "/" .. tostring(c.red) .. ".")
        end
      else
        lua_warn("Unable to resolve the derived potential automatically for couple " .. tostring(c.ox) .. "/" .. tostring(c.red) .. " (insufficient thermodynamic connectivity).")
      end
    end
  end
end

function M.tex_add_template(opts)
  opts = opts or {}
  local raw_items = split_top_level_commas(opts.species or "")
  local element = trim(opts.element or opts.symbol or "")
  if element == "" then element = infer_element_from_species(raw_items) end
  local default_phase = trim(opts.phase_default or "aq")
  if default_phase == "" then default_phase = "aq" end
  local name = trim(opts.name or "")
  local metal = trim(opts.metal or "")
  if metal == "" then metal = element end
  local emit_comments = truthy(opts.comments, true)
  local mode = trim(opts.redox_pairs or "dominant"):lower():gsub("[_-]", " ")
  if mode == "thermodynamic" or mode == "independent" then mode = "dominant" end
  local check_db = truthy(opts.check_database, true)

  template_line("")
  template_line("% ============================================================")
  template_line("% Automatically generated Pourbaix skeleton")
  template_line("% Element : " .. element)
  template_line("% Thermodynamic mode: only independent constants are requested;")
  template_line("% chemical equations associated with pKa values are printed;")
  template_line("% redundant potentials are expressed by thermodynamic cycles.")
  template_line("% The constant RT ln(10)/F is taken as 0.06 V,")
  template_line("% consistently with the plotting engine.")
  template_line("% ============================================================")
  template_line(string.format("\\PourbaixDeclareElement[symbol=%s, name={%s}, metal=%s]", element, name ~= "" and name or element, metal))
  template_line("")

  local species = {}
  for _, raw in ipairs(raw_items) do
    local sp = parse_template_species(raw, element, default_phase)
    if sp.nMetal <= 0 then
      template_line("% WARNING: species `" .. raw .. "` does not appear to contain element " .. element .. ".")
    end
    species[#species+1] = sp
    emit_species_template(element, sp)
  end

  table.sort(species, function(a,b)
    local oa, ob = a.oxState or 0, b.oxState or 0
    if oa == ob then return acid_base_score(a) < acid_base_score(b) end
    return oa < ob
  end)

  local forest = build_acid_base_forest(species)
  local independent_redox = build_independent_redox(species, forest)
  local indep_redox_by_key = {}
  for _,it in ipairs(independent_redox) do indep_redox_by_key[it.key] = it end
  local gexpr = solve_g_expressions(species, forest.edges, independent_redox)

  if truthy(opts.acid_base, true) then
    template_line("")
    template_line("% Independent vertical constants")
    if #forest.edges == 0 and emit_comments then template_line("% No acid/base equilibrium can be inferred automatically.") end
    for _,ab in ipairs(forest.edges) do emit_acidbase_template(element, ab, "TODO") end
  end

  if truthy(opts.redox, true) then
    template_line("")
    template_line("% Redox couples: independent and derived potentials")
    local levels = {}
    local seen = {}
    for _,sp in ipairs(species) do if sp.oxState and not seen[sp.oxState] then levels[#levels+1]=sp.oxState; seen[sp.oxState]=true end end
    table.sort(levels)
    local adjacent = {}; for i=2,#levels do adjacent[levels[i] .. ":" .. levels[i-1]] = true end
    local nR = 0
    for _, ox in ipairs(species) do
      for _, red in ipairs(species) do
        if ox.oxState and red.oxState and ox.oxState > red.oxState then
          local adj = adjacent[ox.oxState .. ":" .. red.oxState] or false
          if should_emit_redox_pair(mode, forest, ox, red, adj) then
            local r = infer_redox(ox, red)
            if r then
              nR = nR + 1
              local key = pair_key(r.ox, r.red)
              local indep = indep_redox_by_key[key]
              local eexpr = nil
              if gexpr[r.ox.id] and gexpr[r.red.id] then
                eexpr = expr_scale(expr_sub(expr_scale(gexpr[r.ox.id], r.nOx or 1), expr_scale(gexpr[r.red.id], r.nRed or 1)), 1/(r.n or 1))
              end
              if indep then
                emit_redox_template(element, r, "TODO", "independent potential: " .. indep.label)
              else
                local formula = eexpr and expr_format(eexpr) or "expression non determinee"
                emit_redox_template(element, r, "auto", "derived potential: " .. redox_label(r) .. " = " .. formula)
              end
            elseif emit_comments then
              template_line("% Unable to balance automatically the couple " .. species_label(ox) .. "/" .. species_label(red) .. ".")
            end
          end
        end
      end
    end
    if nR == 0 and emit_comments then template_line("% No redox couple can be inferred automatically.") end
  end

  template_line("")
  template_line("% Independent numerical data to supply:")
  if truthy(opts.acid_base, true) then
    for _,ab in ipairs(forest.edges) do
      template_line("%   " .. ab.label .. " = TODO")
      template_line("%     equation: " .. acid_base_equation_text(ab))
    end
  end
  if truthy(opts.redox, true) then
    for _,it in ipairs(independent_redox) do template_line("%   " .. it.label .. " = TODO V") end
  end

  if truthy(opts.redox, true) then
    template_line("%")
    template_line("% Redox potentials derived by thermodynamic cycles:")
    local levels = {}; local seen = {}
    for _,sp in ipairs(species) do if sp.oxState and not seen[sp.oxState] then levels[#levels+1]=sp.oxState; seen[sp.oxState]=true end end
    table.sort(levels)
    local adjacent = {}; for i=2,#levels do adjacent[levels[i] .. ":" .. levels[i-1]] = true end
    for _, ox in ipairs(species) do
      for _, red in ipairs(species) do
        if ox.oxState and red.oxState and ox.oxState > red.oxState then
          local adj = adjacent[ox.oxState .. ":" .. red.oxState] or false
          if should_emit_redox_pair(mode, forest, ox, red, adj) then
            local r = infer_redox(ox, red)
            if r then
              local key = pair_key(r.ox, r.red)
              if not indep_redox_by_key[key] then
                local eexpr = nil
                if gexpr[r.ox.id] and gexpr[r.red.id] then
                  eexpr = expr_scale(expr_sub(expr_scale(gexpr[r.ox.id], r.nOx or 1), expr_scale(gexpr[r.red.id], r.nRed or 1)), 1/(r.n or 1))
                end
                template_line("%   " .. redox_label(r) .. " = " .. (eexpr and expr_format(eexpr) or "auto"))
              end
            end
          end
        end
      end
    end
  end

  if check_db and M.data[element] then
    template_line("%")
    template_line("% Check against pourbaix-data.lua:")
    local values = {}
    for _,ab in ipairs(forest.edges) do
      local v = find_database_pKa(element, ab.acid.id, ab.base.id)
      if v ~= nil then values[ab.label] = v end
    end
    for _,it in ipairs(independent_redox) do
      local v = find_database_redox_E0(element, it.r.ox.id, it.r.red.id)
      if v ~= nil then values[it.label] = v end
    end
    local levels = {}; local seen = {}
    for _,sp in ipairs(species) do if sp.oxState and not seen[sp.oxState] then levels[#levels+1]=sp.oxState; seen[sp.oxState]=true end end
    table.sort(levels)
    local adjacent = {}; for i=2,#levels do adjacent[levels[i] .. ":" .. levels[i-1]] = true end
    for _, ox in ipairs(species) do
      for _, red in ipairs(species) do
        if ox.oxState and red.oxState and ox.oxState > red.oxState then
          local adj = adjacent[ox.oxState .. ":" .. red.oxState] or false
          if should_emit_redox_pair(mode, forest, ox, red, adj) then
            local r = infer_redox(ox, red)
            if r and gexpr[r.ox.id] and gexpr[r.red.id] then
              local eexpr = expr_scale(expr_sub(expr_scale(gexpr[r.ox.id], r.nOx or 1), expr_scale(gexpr[r.red.id], r.nRed or 1)), 1/(r.n or 1))
              local calc = expr_eval(eexpr, values)
              local db = find_database_redox_E0(element, r.ox.id, r.red.id)
              if calc and db then
                template_line(string.format("%%   %-28s computed %.6g V ; database %.6g V ; difference %.3g V", redox_label(r), calc, db, calc-db))
              elseif db then
                template_line(string.format("%%   %-28s present in database: %.6g V", redox_label(r), db))
              end
            end
          end
        end
      end
    end
  end

  template_line("% ============================================================")
  template_line("")
end

local function tex_bool(v, default)
  return truthy(v, default)
end

local function get_range_from_opts(elementData, opts)
  local er = elementData.eRange or {-1.0, 2.0}
  return {
    pHmin = number(opts.pHmin, 0),
    pHmax = number(opts.pHmax, 14),
    Emin = number(opts.Emin, er[1]),
    Emax = number(opts.Emax, er[2]),
  }
end

local function sanitize_color(c)
  c = trim(c or "black")
  if c == "" then c = "black" end
  -- xcolor color names and expressions such as black!60 are preserved.
  return c
end

function M.tex_add_pourbaix(opts)
  opts = opts or {}
  local element = trim(opts.element or "Fe")
  local elementData = M.data[element]
  if not elementData then
    tex.error("Unknown Pourbaix element: " .. element, {"Available elements: " .. M.available_elements_string() .. "."})
    return
  end
  local cTrace = number(opts.c, 1e-2)
  if not cTrace or cTrace <= 0 then cTrace = 1e-2 end
  local cConvention = M.normalizeConcentrationConvention(opts.c_convention or opts.concentration_convention or "species")
  local range = get_range_from_opts(elementData, opts)
  local activeSet = M.active_set(elementData, trim(opts.species or "all"))
  local dpH = number(opts.dpH, 0.02)
  local color = sanitize_color(opts.color)
  local lineWidth = trim(opts.line_width or "0.6pt")
  local labelStyle = trim(opts.label_style or "pourbaix label")
  local labelPhase = trim(opts.label_phase or "neutral")
  local domainLabelMode = trim(opts.domain_labels or opts.label_mode or "")
  if domainLabelMode == "" then
    if opts.labels ~= nil and trim(opts.labels) ~= "" then
      domainLabelMode = tex_bool(opts.labels, true) and "formula" or "none"
    else
      domainLabelMode = "formula"
    end
  end
  domainLabelMode = domainLabelMode:lower():gsub("[_-]", " "):gsub("%s+", " ")
  if domainLabelMode == "mhchem" or domainLabelMode == "chemical" or domainLabelMode == "chem" then domainLabelMode = "formula" end
  if domainLabelMode == "letter" or domainLabelMode == "alphabetic" or domainLabelMode == "alphabetical" then domainLabelMode = "letters" end
  if domainLabelMode == "selected domains" or domainLabelMode == "selection" or domainLabelMode == "mixed" or domainLabelMode == "mixte" then domainLabelMode = "selected" end
  local showLabels = not (domainLabelMode == "none" or domainLabelMode == "false" or domainLabelMode == "off" or domainLabelMode == "no")
  local bbox = tex_bool(opts.bbox, true)
  local water = tex_bool(opts.water, false)
  local forget = tex_bool(opts.forget_plot, true)
  local forgetOpt = forget and ", forget plot" or ""

  local frontierLabels = tex_bool(opts.frontier_labels, false)
  local frontierLabelStyle = trim(opts.frontier_label_style or "fill=white, draw, circle, pos=.5, inner sep=1pt, font=\\scriptsize")
  local frontierRedoxLabels = tex_bool(opts.frontier_redox_labels, true)
  local frontierVerticalLabels = tex_bool(opts.frontier_vertical_labels, true)
  local extraXTicks = tex_bool(opts.extra_x_ticks, false)
  local extraYTicks = tex_bool(opts.extra_y_ticks, false)
  local tickSigFigs = math.floor(number(opts.tick_sig_figs, 4) or 4)
  local registryId = trim(opts.id or opts.name or "")

  local fr = M.computeFrontiers(elementData, cTrace, activeSet, range, dpH, cConvention)
  local byName = M.annotateFrontiers(elementData, cTrace, fr, cConvention)
  M.storeFrontierRegistry(registryId, element, cTrace, fr, byName, cConvention)

  local xTickValues, xTickLabels = {}, {}
  local yTickValues, yTickLabels = {}, {}

  if bbox then
    tex.print(string.format("\\addplot[draw=none, mark=none%s] coordinates {(%s,%s) (%s,%s)};",
      forgetOpt, fmt(range.pHmin), fmt(range.Emin), fmt(range.pHmax), fmt(range.Emax)))
  end

  local plotOptions = string.format("draw=%s, line width=%s, mark=none%s", color, lineWidth, forgetOpt)

  for _, s in ipairs(fr.redoxFrontiers) do
    local p1 = s.points[1]
    local p2 = s.points[#s.points]
    if p1 and p2 then
      local node = ""
      if frontierLabels and frontierRedoxLabels then
        node = string.format(" node[%s] {%s}", frontierLabelStyle, s.name)
      end
      tex.print(string.format("\\addplot[%s] coordinates {(%s,%s) (%s,%s)}%s;",
        plotOptions, fmt(p1[1]), fmt(p1[2]), fmt(p2[1]), fmt(p2[2]), node))
      if extraYTicks and s.horizontal and s.intercept then
        append_unique(yTickValues, yTickLabels, s.intercept, fmt_sig(s.intercept, tickSigFigs), 1e-8)
      end
    end
  end

  for _, ab in ipairs(fr.acidBaseFrontiers) do
    local y1 = math.max(range.Emin, ab.eMin)
    local y2 = math.min(range.Emax, ab.eMax)
    if y2 > y1 then
      local node = ""
      if frontierLabels and frontierVerticalLabels then
        node = string.format(" node[%s] {%s}", frontierLabelStyle, ab.name)
      end
      tex.print(string.format("\\addplot[%s] coordinates {(%s,%s) (%s,%s)}%s;",
        plotOptions, fmt(ab.pH), fmt(y1), fmt(ab.pH), fmt(y2), node))
      if extraXTicks then
        append_unique(xTickValues, xTickLabels, ab.pH, fmt_sig(ab.pH, tickSigFigs), 1e-8)
      end
    end
  end

  emit_current_axis_extra_ticks(xTickValues, xTickLabels, yTickValues, yTickLabels, opts.extra_x_tick_style, opts.extra_y_tick_style)

  if water then
    local preset = trim(opts.water_preset or "aere")
    local pO2, pH2 = 0.2, 1
    if preset == "standard" then pO2, pH2 = 1, 1 end
    if preset == "desaere" or preset == "des" then pO2, pH2 = 1e-6, 1 end
    local wf = M.waterFrontiers(pO2, pH2)
    local wstyle = trim(opts.water_style or "draw=black!45, densely dashed, line width=0.45pt, mark=none")
    tex.print(string.format("\\addplot[%s%s] coordinates {(%s,%s) (%s,%s)};",
      wstyle, forgetOpt, fmt(range.pHmin), fmt(wf.H2(range.pHmin)), fmt(range.pHmax), fmt(wf.H2(range.pHmax))))
    tex.print(string.format("\\addplot[%s%s] coordinates {(%s,%s) (%s,%s)};",
      wstyle, forgetOpt, fmt(range.pHmin), fmt(wf.O2(range.pHmin)), fmt(range.pHmax), fmt(wf.O2(range.pHmax))))

    if tex_bool(opts.water_labels, false) then
      local wlabelStyle = trim(opts.water_label_style or "font=\\scriptsize,fill=white,inner sep=1pt")
      local wlabelOffset = number(opts.water_label_offset, 0.10)
      local span = range.pHmax - range.pHmin
      local xLow = number(opts.water_label_pH_low, range.pHmin + 0.14 * span)
      local xHigh = number(opts.water_label_pH_high, range.pHmax - 0.14 * span)
      if xLow < range.pHmin then xLow = range.pHmin end
      if xLow > range.pHmax then xLow = range.pHmax end
      if xHigh < range.pHmin then xHigh = range.pHmin end
      if xHigh > range.pHmax then xHigh = range.pHmax end
      local yH = wf.H2(xLow)
      local yO = wf.O2(xHigh)
      local lowerOx = trim(opts.water_lower_oxidant_label or "H+")
      local lowerRed = trim(opts.water_lower_reductant_label or "H2(aq)")
      local upperOx = trim(opts.water_upper_oxidant_label or "O2(aq)")
      local upperRed = trim(opts.water_upper_reductant_label or "H2O")
      local function water_label_text(s)
        s = trim(s or "")
        if s == "" then return "{}" end
        -- Users normally provide raw mhchem formulas such as O2(g), H2(g), or HO-.
        -- If they provide explicit TeX code, keep it untouched.
        if s:match("\\") then return s end
        return "\\ce{" .. s .. "}"
      end
      lowerOx = water_label_text(lowerOx)
      lowerRed = water_label_text(lowerRed)
      upperOx = water_label_text(upperOx)
      upperRed = water_label_text(upperRed)
      -- Oxidants are placed above the corresponding water line; reductants are
      -- placed below it.  Anchors keep the labels visually separated from the
      -- boundary even when the offset is small.
      tex.print(string.format("\\node[%s,anchor=south] at (axis cs:%s,%s) {%s};", wlabelStyle, fmt(xLow), fmt(yH + wlabelOffset), lowerOx))
      tex.print(string.format("\\node[%s,anchor=north] at (axis cs:%s,%s) {%s};", wlabelStyle, fmt(xLow), fmt(yH - wlabelOffset), lowerRed))
      tex.print(string.format("\\node[%s,anchor=south] at (axis cs:%s,%s) {%s};", wlabelStyle, fmt(xHigh), fmt(yO + wlabelOffset), upperOx))
      tex.print(string.format("\\node[%s,anchor=north] at (axis cs:%s,%s) {%s};", wlabelStyle, fmt(xHigh), fmt(yO - wlabelOffset), upperRed))
    end
  end

  if showLabels then
    local nx = math.floor(number(opts.label_nx, 64))
    local ny = math.floor(number(opts.label_ny, 64))
    local labels = M.label_positions(elementData, cTrace, activeSet, range, nx, ny, nil, cConvention)
    local selectedSet = {}
    local selectedPhase = opts.selected_domain_phase
    if selectedPhase == nil or trim(selectedPhase) == "" then selectedPhase = labelPhase end
    if domainLabelMode == "selected" then
      selectedSet = selected_domain_set(elementData, element, opts.selected_domains or "")
    end
    -- Alphabetic labels use reading order: left to right, then top to bottom
    -- when centroids have the same abscissa.
    table.sort(labels, function(a,b)
      if math.abs((a.pH or 0) - (b.pH or 0)) < 1e-9 then return (a.E or 0) > (b.E or 0) end
      return (a.pH or 0) < (b.pH or 0)
    end)
    local letterIndex = 0
    for i, lab in ipairs(labels) do
      local txt
      if domainLabelMode == "letters" then
        txt = letter_name(i):upper()
      elseif domainLabelMode == "selected" then
        if selectedSet[lab.id] then
          txt = M.label_text(elementData, lab.id, selectedPhase)
        else
          letterIndex = letterIndex + 1
          txt = letter_name(letterIndex):upper()
        end
      else
        txt = M.label_text(elementData, lab.id, labelPhase)
      end
      tex.print(string.format("\\node[%s, text=%s] at (axis cs:%s,%s) {%s};",
        labelStyle, color, fmt(lab.pH), fmt(lab.E), txt))
    end
  end
end

return M
