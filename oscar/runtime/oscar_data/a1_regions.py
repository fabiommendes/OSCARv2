import numpy as np

from oscar.config import mod_regionI, mod_regionJ, mod_sector, ind_final, \
    mod_kindFF, mod_kindLUC, mod_kindGHG, \
    mod_kindCHI, mod_kindAER, mod_kindRF, mod_kindGE, mod_biomeSHR, mod_biomeURB, \
    mod_biomeV3, dty

# ============
# A.1. Regions
# ============
from oscar.data import load_data

nb_regionI = nb_regionJ = 0
mod_regions = {"I": mod_regionI, "J": mod_regionJ}

for X in ["I", "J"]:
    mod_region = mod_regions[X]

    if mod_region == "SRES4":
        region = ["OECD90", "REF", "ASIA", "ALM"]
        region_name = region
        region_color = ["#0000FF", "#660099", "#006600", "#FF6600"]
    elif mod_region == "SRES11":
        region = ["NAM", "WEU", "PAO", "EEU", "FSU", "CPA", "SAS", "PAS", "MEA", "LAM",
                  "AFR"]
        region_name = [
            "North America",
            "Western Europe",
            "Pacific OECD",
            "Central and Eastern Europe",
            "Former Soviet Union",
            "Centrally Planned Asia",
            "South Asia",
            "Pacific Asia",
            "Middle East and North Africa",
            "Latin America and the Caribbean",
            "Sub-Saharan Africa", ]
        region_color = [
            "#000099",
            "#0000FF",
            "#00FFFF",
            "#660099",
            "#FF0099",
            "#006600",
            "#00FF00",
            "#33FF33",
            "#FFFF00",
            "#FF6600",
            "#FF0000", ]
    elif mod_region == "RECCAP*":
        region = ["L1", "L2", "L3", "L4", "L5", "L6", "L7", "L8", "L9", "L10", "LX"]
        region_name = [
            "Africa",
            "Arctic",
            "Oceania",
            "China",
            "S.E. Asia",
            "S. Asia",
            "Europe",
            "N. America",
            "Russia",
            "S. America",
            "Rest", ]
        region_color = ["", "", "", "", "", "", "", "", "", "", ""]
    elif mod_region == "Raupach*":
        region = ["USA", "EU", "Japan", "D1", "FSU", "China", "India", "D2", "D3"]
        region_name = region
        region_color = [
            "#FF0000",
            "#FF6600",
            "#FFFF00",
            "#00FF00",
            "#00FFFF",
            "#0000FF",
            "#FF0099",
            "#660066",
            "#000000", ]
    elif mod_region == "Houghton":
        region = [
            "N. Am.",
            "S. & C. Am.",
            "Europe",
            "N. Afr. & M. East",
            "Trop. Afr.",
            "FSU",
            "China",
            "S. & S.E. Asia",
            "Pacific Dvp.", ]
        region_name = [
            "North America",
            "South & Central America",
            "Europe",
            "North Africa & Middle East",
            "Tropical Africa",
            "Former Soviet Union",
            "China region",
            "South & South-East Asia",
            "Pacific Developed region", ]
        region_color = [
            "#FF0000",
            "#FFFF33",
            "#00FF00",
            "#00FFFF",
            "#0000FF",
            "#FF00FF",
            "#009999",
            "#FF6600",
            "#990066", ]
    elif mod_region == "IMACLIM":
        region = [
            "USA",
            "Canada",
            "Europe",
            "OECD Pacific",
            "CEI",
            "China",
            "India",
            "Brazil",
            "Middle East",
            "Africa",
            "Rest of Asia",
            "Rest of LAM", ]
        region_name = region
        region_color = ["", "", "", "", "", "", "", "", "", "", "", ""]
    elif mod_region == "Kyoto":
        region = ["B", "nB"]
        region_name = ["Annex B", "non- Annex B"]
        region_color = ["#006600", "#CCCC99"]
    elif mod_region == "RCP5":
        region = ["ASIA", "LAM", "MAF", "OECD90", "REF"]
        region_name = [
            "Asia region",
            "Latin America",
            "Middle-East & Africa",
            "OECD countries in 1990",
            "Reforming countries", ]
        region_color = ["", "", "", "", ""]
    elif mod_region == "RCP10*":
        region = [
            "China +",
            "India +",
            "Rest of Asia",
            "Latin America",
            "Middle East",
            "Africa",
            "Western Europe",
            "Northern America",
            "Pacific OECD",
            "Reforming Economies", ]
        region_name = region
        region_color = ["", "", "", "", "", "", "", "", "", ""]
    else:
        region = []
        region_name = []
        region_color = []

    region_index = {}
    path = "data/Regions_GTAP/#DATA.Regions_GTAP.csv"
    TMP = np.array(load_data(path, dtype=object, use='csv')[:,2:])
    for n in range(1, len(TMP)):
        if mod_region in list(TMP[0]):
            region_index[int(TMP[n, 0])] = int(TMP[n, list(TMP[0]).index(mod_region)])
        else:
            region_index[int(TMP[n, 0])] = 0

    if X == 'I':
        regionI = ["n/a"] + region
        regionI_name = ["n/a"] + region_name
        regionI_color = ["0.5"] + region_color
        regionI_index = dict([(0, 0)] + list(region_index.items()))
        nb_regionI = len(regionI)
    else:
        regionJ = ["n/a"] + region
        regionJ_name = ["n/a"] + region_name
        regionJ_color = ["0.5"] + region_color
        regionJ_index = dict([(0, 0)] + list(region_index.items()))
        nb_regionJ = len(regionJ)

    del region, region_name, region_color, region_index

# ============
# A.2. Sectors
# ============


if (mod_sector == "Time") & (300 < ind_final <= 310):
    sector = (
            ["<1850", "1850-1899", "1900-1909", "1910-1919", "1920-1929", "1930-1939",
             "1940-1949", "1950-1959"]
            + [str(t) + "-" + str(t + 4) for t in range(1960, 1990, 5)]
            + [str(t) + "-" + str(t + 1) for t in range(1990, 2000, 2)]
            + [str(t) for t in range(2000, 1700 + ind_final + 1)]
    )
    sector_name = sector
    sector_color = ["0.5" for n in range(len(sector))]
elif (mod_sector == "TimeRCP") & (ind_final == 400):
    sector = ["<2011"] + [str(t + 1) + "-" + str(t + 5) for t in range(2010, 2100, 5)]
    sector_name = sector
    sector_color = ["0.5" for n in range(len(sector))]
else:
    sector = ["n/a"]
    sector_name = sector
    sector_color = ["0.5"]

nb_sector = len(sector)

# ==========
# A.3. Kinds
# ==========

kind = ["n/a"]
kind_name = ["n/a"]
kind_color = ["0.5"]

# fossil CO2
if mod_kindFF == "one":
    kind += ["FF"]
    kind_name += ["Fossil Fuel"]
    kind_color += ["#FF0000"]  # or ['#666666']
    kFF = 1
    kLUC = kFF + 1
elif mod_kindFF == "CDIAC":
    kind += ["FF Solids", "FF Liquids", "FF Gas", "FF Cement", "FF Flaring"]
    kind_name += ["FF Solids", "FF Liquids", "FF Gas", "FF Cement", "FF Flaring"]
    kind_color += ["", "", "", "", ""]
    kFF = 1
    kLUC = kFF + 5
else:
    kFF = 0
    kLUC = 1  # max(0) + 1?

if mod_kindFF == "CDIAC":
    kindFF_index = {"sol": kFF, "liq": kFF + 1, "gas": kFF + 2, "cem": kFF + 3,
                    "fla": kFF + 4}
else:
    kindFF_index = {"sol": kFF, "liq": kFF, "gas": kFF, "cem": kFF, "fla": kFF}

# land-use
if mod_kindLUC == "one":
    kind += ["LULCC"]
    kind_name += ["Land-Use and Land-Cover Change"]
    kind_color += ["#993300"]
    kGHG = kLUC + 1
elif mod_kindLUC == "all":
    kind += ["LUC-CO2", "LUC-BB"]
    kind_name += ["LUC CO2 only", "LUC BB non-CO2"]
    kind_color += []
    kGHG = kLUC + 2
else:
    kLUC = 0
    kGHG = kFF + 1  # max(kFF) + 1?

if mod_kindLUC == "all":
    kindLUC_index = {"CO2": kLUC, "BB": kLUC + 1}
else:
    kindLUC_index = {"CO2": kLUC, "BB": kLUC}

# other GHG
if mod_kindGHG == "one":
    kind += ["non-CO2"]
    kind_name += ["non-CO2"]
    kind_color += ["#FFCC00"]
    kCHI = kGHG + 1
elif mod_kindGHG == "RCP":
    kind += ["CH4", "N2O", "HaloC"]
    kind_name += ["Methane", "Nitrous Oxide", "Halocarbons"]
    kind_color += ["#FF6600", "#FFCC00", "#FF9999"]
    kCHI = kGHG + 3
else:
    kGHG = 0
    kCHI = max(kFF, kLUC) + 1

if mod_kindGHG == "RCP":
    kindGHG_index = {"CH4": kGHG, "N2O": kGHG + 1, "HFC": kGHG + 2, "PFC": kGHG + 2,
                     "ODS": kGHG + 2}
else:
    kindGHG_index = {"CH4": kGHG, "N2O": kGHG, "HFC": kGHG, "PFC": kGHG, "ODS": kGHG}

# active species
if mod_kindCHI == "one":
    kind += ["OzPrec."]
    kind_name += ["Ozone Precursors"]
    kind_color += ["#66FF66"]
    kAER = kCHI + 1
elif mod_kindCHI == "all":
    kind += ["NOx", "CO", "NMVOC"]
    kind_name += ["Nitrogen Oxides", "Carbon Monoxide",
                  "Non-Methane Volatile Organic Compounds"]
    kind_color += ["#66FF66", "#009999", "#006600"]
    kAER = kCHI + 3
else:
    kCHI = 0
    kAER = max(kFF, kLUC, kGHG) + 1

if mod_kindCHI == "all":
    kindCHI_index = {"NOX": kCHI, "CO": kCHI + 1, "VOC": kCHI + 2}
else:
    kindCHI_index = {"NOX": kCHI, "CO": kCHI, "VOC": kCHI}

# aerosols
if mod_kindAER == "one":
    kind += ["AER"]
    kind_name += ["Aerosols"]
    kind_color += ["#0000FF"]
    kRF = kAER + 1
elif mod_kindAER == "all":
    kind += ["SO2", "NH3", "OC", "BC"]
    kind_name += ["Sulfur Dioxide", "Ammonia", "Organic Carbon", "Black Carbon"]
    kind_color += ["#00CCFF", "#0000FF", "#660099", "#CC0066"]
    kRF = kAER + 4
else:
    kAER = 0
    kRF = max(kFF, kLUC, kGHG, kCHI) + 1

if mod_kindAER == "all":
    kindAER_index = {"SO2": kAER, "NH3": kAER + 1, "OC": kAER + 2, "BC": kAER + 3}
else:
    kindAER_index = {"SO2": kAER, "NH3": kAER, "OC": kAER, "BC": kAER}

# other radiative forcings
if mod_kindRF == "one":
    kind += ["RFother"]
    kind_name += ["Other RF"]
    kind_color += ["#999999"]
    kGE = kRF + 1
elif mod_kindRF == "two":
    kind += ["RFant", "RFnat"]
    kind_name += ["Anthropogenic RF", "Natural RF"]
    kind_color += ["", ""]
    kGE = kRF + 2
elif mod_kindRF == "all":
    kind += ["RFcon", "RFsol", "RFvol"]
    kind_name += ["Contrails RF", "Solar RF", "Volcanoes RF"]
    kind_color += ["", "", ""]
    kGE = kRF + 3
else:
    kRF = 0
    kGE = max(kFF, kLUC, kGHG, kCHI, kAER) + 1

if mod_kindRF == "all":
    kindRF_index = {"RFcon": kRF, "RFsol": kRF + 1, "RFvol": kRF + 2}
elif mod_kindRF == "two":
    kindRF_index = {"RFcon": kRF, "RFsol": kRF + 1, "RFvol": kRF + 1}
else:
    kindRF_index = {"RFcon": kRF, "RFsol": kRF, "RFvol": kRF}

# Geoengineering
if mod_kindGE == "PUP":
    kind += ["AFO", "CCS", "ALB", "AER"]
    kind_name += ["Aforestation", "Carbon Capture and Storage", "Surface Albedo",
                  "Sulfate Aerosols"]
    kind_color += ["", "", "", ""]
else:
    kGE = 0

nb_kind = len(kind)

# ===========
# A.4. Biomes
# ===========

if (mod_biomeSHR == "w/FOR") & (mod_biomeURB == "w/DES"):
    biome = ["DES+", "FOR+", "GRA", "CRO", "PAS"]
    biome_name = ["Desert & Urban", "Forest & Shrubland", "Grassland", "Cropland",
                  "Pasture"]
    biome_color = ["", "", "", "", ""]
    biome_index = {"des": 0, "for": 1, "shr": 1, "gra": 2, "cro": 3, "pas": 4, "urb": 0}
elif (mod_biomeSHR == "w/FOR") & (mod_biomeURB == "URB"):
    biome = ["DES", "FOR+", "GRA", "CRO", "PAS", "URB"]
    biome_name = ["Desert", "Forest & Shrubland", "Grassland", "Cropland", "Pasture",
                  "Urban"]
    biome_color = ["", "", "", "", "", ""]
    biome_index = {"des": 0, "for": 1, "shr": 1, "gra": 2, "cro": 3, "pas": 4, "urb": 5}
elif (mod_biomeSHR == "w/GRA") & (mod_biomeURB == "w/DES"):
    biome = ["DES+", "FOR", "GRA+", "CRO", "PAS"]
    biome_name = ["Desert & Urban", "Forest", "Grassland & Shrubland", "Cropland",
                  "Pasture"]
    biome_color = ["", "", "", "", ""]
    biome_index = {"des": 0, "for": 1, "shr": 2, "gra": 2, "cro": 3, "pas": 4, "urb": 0}
elif (mod_biomeSHR == "w/GRA") & (mod_biomeURB == "URB"):
    biome = ["DES", "FOR", "GRA+", "CRO", "PAS", "URB"]
    biome_name = ["Desert", "Forest", "Grassland & Shrubland", "Cropland", "Pasture",
                  "Urban"]
    biome_color = ["", "", "", "", "", ""]
    biome_index = {"des": 0, "for": 1, "shr": 2, "gra": 2, "cro": 3, "pas": 4, "urb": 5}
elif (mod_biomeSHR == "SHR") & (mod_biomeURB == "w/DES"):
    biome = ["DES+", "FOR", "SHR", "GRA", "CRO", "PAS"]
    biome_name = ["Desert & Urban", "Forest", "Shrubland", "Grassland", "Cropland",
                  "Pasture"]
    biome_color = ["", "", "", "", "", ""]
    biome_index = {"des": 0, "for": 1, "shr": 2, "gra": 3, "cro": 4, "pas": 5, "urb": 0}
elif (mod_biomeSHR == "SHR") & (mod_biomeURB == "URB"):
    biome = ["DES", "FOR", "SHR", "GRA", "CRO", "PAS", "URB"]
    biome_name = ["Desert", "Forest", "Shrubland", "Grassland", "Cropland", "Pasture",
                  "Urban"]
    biome_color = ["", "", "", "", "", "", ""]
    biome_index = {"des": 0, "for": 1, "shr": 2, "gra": 3, "cro": 4, "pas": 5, "urb": 6}
else:
    biome = ["all"]
    biome_name = biome
    biome_color = []
    biome_index = {"des": 0, "for": 0, "shr": 0, "gra": 0, "cro": 0, "pas": 0, "urb": 0}

# for OSCAR v3 parameters
if mod_biomeV3:
    biome = ["FOR", "GRA+", "CRO", "PAS", "URB"]
    biome_name = ["Forest", "Non-Forest", "Cropland", "Pasture", "Urban"]
    biome_color = ["", "", "", "", ""]
    biome_index = {"des": 1, "for": 0, "shr": 1, "gra": 1, "cro": 2, "pas": 3, "urb": 4}

nb_biome = len(biome)

# =========
# A.5. Halo
# =========

HFC = [
    "HFC23",
    "HFC32",
    "HFC125",
    "HFC134a",
    "HFC143a",
    "HFC152a",
    "HFC227ea",
    "HFC236fa",
    "HFC245fa",
    "HFC365mfc",
    "HFC4310mee",
]
PFC = ["SF6", "NF3", "CF4", "C2F6", "C3F8", "cC4F8", "C4F10", "C5F12", "C6F14", "C7F16"]
ODS = [
    "CFC11",
    "CFC12",
    "CFC113",
    "CFC114",
    "CFC115",
    "CCl4",
    "CH3CCl3",
    "HCFC22",
    "HCFC141b",
    "HCFC142b",
    "Halon1211",
    "Halon1202",
    "Halon1301",
    "Halon2402",
    "CH3Br",
    "CH3Cl",
]

nb_HFC = len(HFC)
nb_PFC = len(PFC)
nb_ODS = len(ODS)
