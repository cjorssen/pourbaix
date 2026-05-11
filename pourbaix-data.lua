-- This file is part of the pourbaix package.
-- Version: 1.0 (2026-05-10).
-- License: LaTeX Project Public License v1.3c or later.
-- Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>.

-- Built-in thermodynamic data for the Pourbaix engine.
return {
  Fe = {
    name = "Iron",
    metal = "Fe",
    species = {
      Fe = {
        formula = "Fe",
        phase = "s",
        oxState = 0,
        color = "#9a9a9a"
      },
      ["Fe2+"] = {
        formula = "Fe^{2+}",
        phase = "aq",
        oxState = 2,
        color = "#9bcf8d"
      },
      ["Fe(OH)2"] = {
        formula = "Fe(OH)_2",
        phase = "s",
        oxState = 2,
        color = "#3d7a3d"
      },
      ["Fe3+"] = {
        formula = "Fe^{3+}",
        phase = "aq",
        oxState = 3,
        color = "#d4a04a"
      },
      ["Fe(OH)3"] = {
        formula = "Fe(OH)_3",
        phase = "s",
        oxState = 3,
        color = "#b04a2a"
      }
    },
    acidBaseEquilibria = {
      {
        sideAcid = "Fe2+",
        sideBase = "Fe(OH)2",
        nH = 2,
        nH2O = 2,
        pKa = 12.9
      },
      {
        sideAcid = "Fe3+",
        sideBase = "Fe(OH)3",
        nH = 3,
        nH2O = 3,
        pKa = 4
      }
    },
    redoxCouples = {
      {
        ox = "Fe2+",
        red = "Fe",
        E0 = -0.44,
        n = 2,
        nH = 0,
        nH2O = 0
      },
      {
        ox = "Fe(OH)2",
        red = "Fe",
        E0 = -0.053,
        n = 2,
        nH = 2,
        nH2O = 2
      },
      {
        ox = "Fe3+",
        red = "Fe2+",
        E0 = 0.77,
        n = 1,
        nH = 0,
        nH2O = 0
      },
      {
        ox = "Fe(OH)3",
        red = "Fe2+",
        E0 = 1.01,
        n = 1,
        nH = 3,
        nH2O = 3
      },
      {
        ox = "Fe(OH)3",
        red = "Fe(OH)2",
        E0 = 0.236,
        n = 1,
        nH = 1,
        nH2O = 1
      },
      {
        ox = "Fe3+",
        red = "Fe",
        E0 = -0.037,
        n = 3,
        nH = 0,
        nH2O = 0
      },
      {
        ox = "Fe(OH)3",
        red = "Fe",
        E0 = 0.043,
        n = 3,
        nH = 3,
        nH2O = 3
      }
    }
  },
  Cu = {
    name = "Copper",
    metal = "Cu",
    species = {
      Cu = {
        formula = "Cu",
        phase = "s",
        oxState = 0,
        color = "#c87653"
      },
      ["Cu2+"] = {
        formula = "Cu^{2+}",
        phase = "aq",
        oxState = 2,
        color = "#3a8fc9"
      },
      ["Cu(OH)2"] = {
        formula = "Cu(OH)_2",
        phase = "s",
        oxState = 2,
        color = "#5da7d4"
      },
      Cu2O = {
        formula = "Cu_2O",
        phase = "s",
        oxState = 1,
        color = "#b04020",
        nMetal = 2
      }
    },
    acidBaseEquilibria = {
      {
        sideAcid = "Cu2+",
        sideBase = "Cu(OH)2",
        nH = 2,
        nH2O = 2,
        pKa = 8.3
      }
    },
    redoxCouples = {
      {
        ox = "Cu2+",
        red = "Cu",
        E0 = 0.34,
        n = 2,
        nH = 0,
        nH2O = 0
      },
      {
        ox = "Cu(OH)2",
        red = "Cu",
        E0 = 0.589,
        n = 2,
        nH = 2,
        nH2O = 2
      },
      {
        ox = "Cu2O",
        red = "Cu",
        E0 = 0.46,
        n = 2,
        nH = 2,
        nH2O = 1,
        nOx = 1,
        nRed = 2
      },
      {
        ox = "Cu2+",
        red = "Cu2O",
        E0 = 0.22,
        n = 2,
        nH = 0,
        nH_right = 2,
        nH2O = 0,
        nH2O_left = 1,
        nOx = 2,
        nRed = 1
      },
      {
        ox = "Cu(OH)2",
        red = "Cu2O",
        E0 = 0.718,
        n = 2,
        nH = 2,
        nH2O = 3,
        nOx = 2,
        nRed = 1
      }
    }
  },
  Zn = {
    name = "Zinc",
    metal = "Zn",
    species = {
      Zn = {
        formula = "Zn",
        phase = "s",
        oxState = 0,
        color = "#8a8a98"
      },
      ["Zn2+"] = {
        formula = "Zn^{2+}",
        phase = "aq",
        oxState = 2,
        color = "#5a9fb8"
      },
      ["Zn(OH)2"] = {
        formula = "Zn(OH)_2",
        phase = "s",
        oxState = 2,
        color = "#b3a4cf"
      },
      ["ZnO22-"] = {
        formula = "ZnO_2^{2-}",
        phase = "aq",
        oxState = 2,
        color = "#7fa67a"
      }
    },
    acidBaseEquilibria = {
      {
        sideAcid = "Zn2+",
        sideBase = "Zn(OH)2",
        nH = 2,
        nH2O = 2,
        pKa = 11
      },
      {
        sideAcid = "Zn(OH)2",
        sideBase = "ZnO22-",
        nH = 2,
        nH2O = 0,
        pKa = 28
      }
    },
    redoxCouples = {
      {
        ox = "Zn2+",
        red = "Zn",
        E0 = -0.76,
        n = 2,
        nH = 0,
        nH2O = 0
      },
      {
        ox = "Zn(OH)2",
        red = "Zn",
        E0 = -0.43,
        n = 2,
        nH = 2,
        nH2O = 2
      },
      {
        ox = "ZnO22-",
        red = "Zn",
        E0 = 0.41,
        n = 2,
        nH = 4,
        nH2O = 2
      }
    }
  },
  Al = {
    name = "Aluminum",
    metal = "Al",
    eRange = {
      -2.5,
      1.5
    },
    species = {
      Al = {
        formula = "Al",
        phase = "s",
        oxState = 0,
        color = "#9aa8b8"
      },
      ["Al3+"] = {
        formula = "Al^{3+}",
        phase = "aq",
        oxState = 3,
        color = "#d4c08e"
      },
      ["Al(OH)3"] = {
        formula = "Al(OH)_3",
        phase = "s",
        oxState = 3,
        color = "#bcd0d0"
      },
      ["Al(OH)4-"] = {
        formula = "Al(OH)_4^{-}",
        phase = "aq",
        oxState = 3,
        color = "#c8a4c0"
      }
    },
    acidBaseEquilibria = {
      {
        sideAcid = "Al3+",
        sideBase = "Al(OH)3",
        nH = 3,
        nH2O = 3,
        pKa = 9
      },
      {
        sideAcid = "Al(OH)3",
        sideBase = "Al(OH)4-",
        nH = 1,
        nH2O = 1,
        pKa = 14
      }
    },
    redoxCouples = {
      {
        ox = "Al3+",
        red = "Al",
        E0 = -1.66,
        n = 3,
        nH = 0,
        nH2O = 0
      },
      {
        ox = "Al(OH)3",
        red = "Al",
        E0 = -1.48,
        n = 3,
        nH = 3,
        nH2O = 3
      },
      {
        ox = "Al(OH)4-",
        red = "Al",
        E0 = -1.2,
        n = 3,
        nH = 4,
        nH2O = 4
      }
    }
  },
  Ag = {
    name = "Silver",
    metal = "Ag",
    species = {
      Ag = {
        formula = "Ag",
        phase = "s",
        oxState = 0,
        color = "#b8b8c4"
      },
      ["Ag+"] = {
        formula = "Ag^{+}",
        phase = "aq",
        oxState = 1,
        color = "#d8c890"
      },
      Ag2O = {
        formula = "Ag_2O",
        phase = "s",
        oxState = 1,
        color = "#5a4030",
        nMetal = 2
      }
    },
    acidBaseEquilibria = {
      {
        sideAcid = "Ag+",
        sideBase = "Ag2O",
        nAcid = 2,
        nBase = 1,
        nH = 2,
        nH2O = 1,
        pKa = 12.2
      }
    },
    redoxCouples = {
      {
        ox = "Ag+",
        red = "Ag",
        E0 = 0.8,
        n = 1,
        nH = 0,
        nH2O = 0
      },
      {
        ox = "Ag2O",
        red = "Ag",
        E0 = 1.166,
        n = 2,
        nH = 2,
        nH2O = 1,
        nOx = 1,
        nRed = 2
      }
    }
  },
  Pb = {
    name = "Lead",
    metal = "Pb",
    species = {
      Pb = {
        formula = "Pb",
        phase = "s",
        oxState = 0,
        color = "#5a5a64"
      },
      ["Pb2+"] = {
        formula = "Pb^{2+}",
        phase = "aq",
        oxState = 2,
        color = "#c4ad7a"
      },
      ["Pb(OH)2"] = {
        formula = "Pb(OH)_2",
        phase = "s",
        oxState = 2,
        color = "#d0d4dc"
      },
      ["HPbO2-"] = {
        formula = "HPbO_2^{-}",
        phase = "aq",
        oxState = 2,
        color = "#d8a8a8"
      }
    },
    acidBaseEquilibria = {
      {
        sideAcid = "Pb2+",
        sideBase = "Pb(OH)2",
        nH = 2,
        nH2O = 2,
        pKa = 13
      },
      {
        sideAcid = "Pb(OH)2",
        sideBase = "HPbO2-",
        nH = 1,
        nH2O = 0,
        pKa = 14.5
      }
    },
    redoxCouples = {
      {
        ox = "Pb2+",
        red = "Pb",
        E0 = -0.13,
        n = 2,
        nH = 0,
        nH2O = 0
      },
      {
        ox = "Pb(OH)2",
        red = "Pb",
        E0 = 0.26,
        n = 2,
        nH = 2,
        nH2O = 2
      },
      {
        ox = "HPbO2-",
        red = "Pb",
        E0 = 0.695,
        n = 2,
        nH = 3,
        nH2O = 2
      }
    }
  },
  Ni = {
    name = "Nickel",
    metal = "Ni",
    species = {
      Ni = {
        formula = "Ni",
        phase = "s",
        oxState = 0,
        color = "#9bb09e"
      },
      ["Ni2+"] = {
        formula = "Ni^{2+}",
        phase = "aq",
        oxState = 2,
        color = "#7bb085"
      },
      ["Ni(OH)2"] = {
        formula = "Ni(OH)_2",
        phase = "s",
        oxState = 2,
        color = "#b8d4be"
      }
    },
    acidBaseEquilibria = {
      {
        sideAcid = "Ni2+",
        sideBase = "Ni(OH)2",
        nH = 2,
        nH2O = 2,
        pKa = 12.7
      }
    },
    redoxCouples = {
      {
        ox = "Ni2+",
        red = "Ni",
        E0 = -0.25,
        n = 2,
        nH = 0,
        nH2O = 0
      },
      {
        ox = "Ni(OH)2",
        red = "Ni",
        E0 = 0.131,
        n = 2,
        nH = 2,
        nH2O = 2
      }
    }
  },
  Au = {
    name = "Gold",
    metal = "Au",
    species = {
      Au = {
        formula = "Au",
        phase = "s",
        oxState = 0,
        color = "#e0b441"
      },
      ["Au3+"] = {
        formula = "Au^{3+}",
        phase = "aq",
        oxState = 3,
        color = "#d8a060"
      }
    },
    acidBaseEquilibria = {},
    redoxCouples = {
      {
        ox = "Au3+",
        red = "Au",
        E0 = 1.5,
        n = 3,
        nH = 0,
        nH2O = 0
      }
    }
  },
  Cr = {
    name = "Chromium",
    metal = "Cr",
    species = {
      Cr = {
        formula = "Cr",
        phase = "s",
        oxState = 0,
        color = "#a8a8b8"
      },
      ["Cr3+"] = {
        formula = "Cr^{3+}",
        phase = "aq",
        oxState = 3,
        color = "#5fa074"
      },
      ["Cr(OH)3"] = {
        formula = "Cr(OH)_3",
        phase = "s",
        oxState = 3,
        color = "#7eb594"
      },
      ["CrO4 2-"] = {
        formula = "CrO_4^{2-}",
        phase = "aq",
        oxState = 6,
        color = "#e8d04a"
      },
      ["Cr2O7 2-"] = {
        formula = "Cr_2O_7^{2-}",
        phase = "aq",
        oxState = 6,
        color = "#e88a2c",
        nMetal = 2
      }
    },
    acidBaseEquilibria = {
      {
        sideAcid = "Cr3+",
        sideBase = "Cr(OH)3",
        nH = 3,
        nH2O = 3,
        pKa = 11.8
      },
      {
        sideAcid = "Cr2O7 2-",
        sideBase = "CrO4 2-",
        nAcid = 1,
        nBase = 2,
        nH = 2,
        nH2O = 1,
        pKa = 14.6
      }
    },
    redoxCouples = {
      {
        ox = "Cr3+",
        red = "Cr",
        E0 = -0.74,
        n = 3,
        nH = 0,
        nH2O = 0
      },
      {
        ox = "Cr(OH)3",
        red = "Cr",
        E0 = -0.504,
        n = 3,
        nH = 3,
        nH2O = 3
      },
      {
        ox = "Cr2O7 2-",
        red = "Cr3+",
        E0 = 1.33,
        n = 6,
        nH = 14,
        nH2O = 7,
        nOx = 1,
        nRed = 2
      },
      {
        ox = "Cr2O7 2-",
        red = "Cr(OH)3",
        E0 = 1.094,
        n = 6,
        nH = 8,
        nH2O = 1,
        nOx = 1,
        nRed = 2
      },
      {
        ox = "CrO4 2-",
        red = "Cr3+",
        E0 = 1.476,
        n = 3,
        nH = 8,
        nH2O = 4
      },
      {
        ox = "CrO4 2-",
        red = "Cr(OH)3",
        E0 = 1.24,
        n = 3,
        nH = 5,
        nH2O = 1
      },
      {
        ox = "Cr2O7 2-",
        red = "Cr",
        E0 = 0.295,
        n = 12,
        nH = 14,
        nH2O = 7,
        nOx = 1,
        nRed = 2
      },
      {
        ox = "CrO4 2-",
        red = "Cr",
        E0 = 0.368,
        n = 6,
        nH = 8,
        nH2O = 4
      }
    }
  },
  Mn = {
    name = "Manganese",
    metal = "Mn",
    eRange = {
      -1.5,
      2.4
    },
    species = {
      Mn = {
        formula = "Mn",
        phase = "s",
        oxState = 0,
        color = "#a09898"
      },
      ["Mn2+"] = {
        formula = "Mn^{2+}",
        phase = "aq",
        oxState = 2,
        color = "#e8b8c8"
      },
      ["Mn(OH)2"] = {
        formula = "Mn(OH)_2",
        phase = "s",
        oxState = 2,
        color = "#e8d4d8"
      },
      MnO2 = {
        formula = "MnO_2",
        phase = "s",
        oxState = 4,
        color = "#3a342e"
      },
      ["MnO4 2-"] = {
        formula = "MnO_4^{2-}",
        phase = "aq",
        oxState = 6,
        color = "#3a8050"
      },
      ["MnO4-"] = {
        formula = "MnO_4^{-}",
        phase = "aq",
        oxState = 7,
        color = "#7a3a98"
      }
    },
    acidBaseEquilibria = {
      {
        sideAcid = "Mn2+",
        sideBase = "Mn(OH)2",
        nH = 2,
        nH2O = 2,
        pKa = 15
      }
    },
    redoxCouples = {
      {
        ox = "Mn2+",
        red = "Mn",
        E0 = -1.18,
        n = 2,
        nH = 0,
        nH2O = 0
      },
      {
        ox = "Mn(OH)2",
        red = "Mn",
        E0 = -0.73,
        n = 2,
        nH = 2,
        nH2O = 2
      },
      {
        ox = "MnO2",
        red = "Mn2+",
        E0 = 1.23,
        n = 2,
        nH = 4,
        nH2O = 2
      },
      {
        ox = "MnO2",
        red = "Mn(OH)2",
        E0 = 0.78,
        n = 2,
        nH = 2,
        nH2O = 2
      },
      {
        ox = "MnO4 2-",
        red = "MnO2",
        E0 = 2.26,
        n = 2,
        nH = 4,
        nH2O = 2
      },
      {
        ox = "MnO4-",
        red = "MnO4 2-",
        E0 = 0.56,
        n = 1,
        nH = 0,
        nH2O = 0
      },
      {
        ox = "MnO4-",
        red = "MnO2",
        E0 = 1.7,
        n = 3,
        nH = 4,
        nH2O = 2
      },
      {
        ox = "MnO4-",
        red = "Mn2+",
        E0 = 1.51,
        n = 5,
        nH = 8,
        nH2O = 4
      },
      {
        ox = "MnO4-",
        red = "Mn(OH)2",
        E0 = 1.329,
        n = 5,
        nH = 6,
        nH2O = 4
      },
      {
        ox = "MnO2",
        red = "Mn",
        E0 = 0.025,
        n = 4,
        nH = 4,
        nH2O = 2
      },
      {
        ox = "MnO4-",
        red = "Mn",
        E0 = 0.741,
        n = 7,
        nH = 8,
        nH2O = 4
      },
      {
        ox = "MnO4 2-",
        red = "Mn",
        E0 = 0.771,
        n = 6,
        nH = 8,
        nH2O = 4
      }
    }
  },
  Cl = {
    name = "Chlorine",
    metal = "Cl-",
    isMetal = false,
    species = {
      ["Cl-"] = {
        formula = "Cl^{-}",
        phase = "aq",
        oxState = -1,
        color = "#d8e8b8"
      },
      Cl2 = {
        formula = "Cl_2",
        phase = "g",
        oxState = 0,
        color = "#bcd866",
        nMetal = 2
      },
      HClO = {
        formula = "HClO",
        phase = "aq",
        oxState = 1,
        color = "#f0e88c"
      },
      ["ClO-"] = {
        formula = "ClO^{-}",
        phase = "aq",
        oxState = 1,
        color = "#e8d870"
      }
    },
    acidBaseEquilibria = {
      {
        sideAcid = "HClO",
        sideBase = "ClO-",
        nH = 1,
        nH2O = 0,
        pKa = 7.5
      }
    },
    redoxCouples = {
      {
        ox = "Cl2",
        red = "Cl-",
        E0 = 1.36,
        n = 2,
        nH = 0,
        nH2O = 0,
        nOx = 1,
        nRed = 2
      },
      {
        ox = "HClO",
        red = "Cl2",
        E0 = 1.63,
        n = 2,
        nH = 2,
        nH2O = 2,
        nOx = 2,
        nRed = 1
      },
      {
        ox = "ClO-",
        red = "Cl2",
        E0 = 2.08,
        n = 2,
        nH = 4,
        nH2O = 2,
        nOx = 2,
        nRed = 1
      },
      {
        ox = "HClO",
        red = "Cl-",
        E0 = 1.495,
        n = 2,
        nH = 1,
        nH2O = 1
      },
      {
        ox = "ClO-",
        red = "Cl-",
        E0 = 1.72,
        n = 2,
        nH = 2,
        nH2O = 1
      }
    }
  }
}
