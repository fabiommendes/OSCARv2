import os
import csv

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

from scipy.optimize import fmin,fsolve
from scipy.special import gammainc


################################################################################
# OSCAR LOADD FILE
################################################################################

print 'LOADING: DRIVERS'

##################################################
#   A. VECTORS
##################################################

# ============
# A.1. Regions
# ============

for X in ['I', 'J']:

    exec('mod_region = mod_region' + X)

    if (mod_region == 'SRES4'):
        region = ['OECD90', 'REF', 'ASIA', 'ALM']
        region_name = region
        region_color = ['#0000FF', '#660099', '#006600', '#FF6600']
    elif (mod_region == 'SRES11'):
        region = ['NAM', 'WEU', 'PAO', 'EEU', 'FSU', 'CPA', 'SAS', 'PAS', 'MEA', 'LAM',
                  'AFR']
        region_name = ['North America', 'Western Europe', 'Pacific OECD',
                       'Central and Eastern Europe', 'Former Soviet Union',
                       'Centrally Planned Asia', 'South Asia', 'Pacific Asia',
                       'Middle East and North Africa', 'Latin America and the Caribbean',
                       'Sub-Saharan Africa']
        region_color = ['#000099', '#0000FF', '#00FFFF', '#660099', '#FF0099', '#006600',
                        '#00FF00', '#33FF33', '#FFFF00', '#FF6600', '#FF0000']
    elif (mod_region == 'RECCAP*'):
        region = ['L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'L8', 'L9', 'L10', 'LX']
        region_name = ['Africa', 'Arctic', 'Oceania', 'China', 'S.E. Asia', 'S. Asia',
                       'Europe', 'N. America', 'Russia', 'S. America', 'Rest']
        region_color = ['', '', '', '', '', '', '', '', '', '', '']
    elif (mod_region == 'Raupach*'):
        region = ['USA', 'EU', 'Japan', 'D1', 'FSU', 'China', 'India', 'D2', 'D3']
        region_name = region
        region_color = ['#FF0000', '#FF6600', '#FFFF00', '#00FF00', '#00FFFF', '#0000FF',
                        '#FF0099', '#660066', '#000000']
    elif (mod_region == 'Houghton'):
        region = ['N. Am.', 'S. & C. Am.', 'Europe', 'N. Afr. & M. East', 'Trop. Afr.',
                  'FSU', 'China', 'S. & S.E. Asia', 'Pacific Dvp.']
        region_name = ['North America', 'South & Central America', 'Europe',
                       'North Africa & Middle East', 'Tropical Africa',
                       'Former Soviet Union', 'China region', 'South & South-East Asia',
                       'Pacific Developed region']
        region_color = ['#FF0000', '#FFFF33', '#00FF00', '#00FFFF', '#0000FF', '#FF00FF',
                        '#009999', '#FF6600', '#990066']
    elif (mod_region == 'IMACLIM'):
        region = ['USA', 'Canada', 'Europe', 'OECD Pacific', 'CEI', 'China', 'India',
                  'Brazil', 'Middle East', 'Africa', 'Rest of Asia', 'Rest of LAM']
        region_name = region
        region_color = ['', '', '', '', '', '', '', '', '', '', '', '']
    elif (mod_region == 'Kyoto'):
        region = ['B', 'nB']
        region_name = ['Annex B', 'non- Annex B']
        region_color = ['#006600', '#CCCC99']
    elif (mod_region == 'RCP5'):
        region = ['ASIA', 'LAM', 'MAF', 'OECD90', 'REF']
        region_name = ['Asia region', 'Latin America', 'Middle-East & Africa',
                       'OECD countries in 1990', 'Reforming countries']
        region_color = ['', '', '', '', '']
    elif (mod_region == 'RCP10*'):
        region = ['China +', 'India +', 'Rest of Asia', 'Latin America', 'Middle East',
                  'Africa', 'Western Europe', 'Northern America', 'Pacific OECD',
                  'Reforming Economies']
        region_name = region
        region_color = ['', '', '', '', '', '', '', '', '', '']
    else:
        region = []
        region_name = []
        region_color = []

    region_index = {}
    TMP = np.array([line for line in
                    csv.reader(open('data/Regions_GTAP/#DATA.Regions_GTAP.csv', 'r'))])[:,
          2:]
    for n in range(1, len(TMP)):
        if mod_region in list(TMP[0]):
            region_index[int(TMP[n, 0])] = int(TMP[n, list(TMP[0]).index(mod_region)])
        else:
            region_index[int(TMP[n, 0])] = 0

    exec('region' + X + ' = ["n/a"]+region')
    exec('region' + X + '_name = ["n/a"]+region_name')
    exec('region' + X + '_color = ["0.5"]+region_color')
    exec('region' + X + '_index = dict([(0,0)]+region_index.items())')
    exec('nb_region' + X + ' = len(region' + X + ')')

    del region, region_name, region_color, region_index

# ============
# A.2. Sectors
# ============

if (mod_sector == 'Time') & (300 < ind_final <= 310):
    sector = ['<1850', '1850-1899', '1900-1909', '1910-1919', '1920-1929', '1930-1939',
              '1940-1949', '1950-1959'] \
             + [str(t) + '-' + str(t + 4) for t in range(1960, 1990, 5)] \
             + [str(t) + '-' + str(t + 1) for t in range(1990, 2000, 2)] \
             + [str(t) for t in range(2000, 1700 + ind_final + 1)]
    sector_name = sector
    sector_color = ['0.5' for n in range(len(sector))]
elif (mod_sector == 'TimeRCP') & (ind_final == 400):
    sector = ['<2011'] + [str(t + 1) + '-' + str(t + 5) for t in range(2010, 2100, 5)]
    sector_name = sector
    sector_color = ['0.5' for n in range(len(sector))]
else:
    sector = ['n/a']
    sector_name = sector
    sector_color = ['0.5']

nb_sector = len(sector)

# ==========
# A.3. Kinds
# ==========

kind = ['n/a']
kind_name = ['n/a']
kind_color = ['0.5']

# fossil CO2
if (mod_kindFF == 'one'):
    kind += ['FF']
    kind_name += ['Fossil Fuel']
    kind_color += ['#FF0000']  # or ['#666666']
    kFF = 1
    kLUC = kFF + 1
elif (mod_kindFF == 'CDIAC'):
    kind += ['FF Solids', 'FF Liquids', 'FF Gas', 'FF Cement', 'FF Flaring']
    kind_name += ['FF Solids', 'FF Liquids', 'FF Gas', 'FF Cement', 'FF Flaring']
    kind_color += ['', '', '', '', '']
    kFF = 1
    kLUC = kFF + 5
else:
    kFF = 0
    kLUC = max(0) + 1

if (mod_kindFF == 'CDIAC'):
    kindFF_index = {'sol': kFF, 'liq': kFF + 1, 'gas': kFF + 2, 'cem': kFF + 3,
                    'fla': kFF + 4}
else:
    kindFF_index = {'sol': kFF, 'liq': kFF, 'gas': kFF, 'cem': kFF, 'fla': kFF}

# land-use
if (mod_kindLUC == 'one'):
    kind += ['LULCC']
    kind_name += ['Land-Use and Land-Cover Change']
    kind_color += ['#993300']
    kGHG = kLUC + 1
elif (mod_kindLUC == 'all'):
    kind += ['LUC-CO2', 'LUC-BB']
    kind_name += ['LUC CO2 only', 'LUC BB non-CO2']
    kind_color += []
    kGHG = kLUC + 2
else:
    kLUC = 0
    kGHG = max(kFF) + 1

if (mod_kindLUC == 'all'):
    kindLUC_index = {'CO2': kLUC, 'BB': kLUC + 1}
else:
    kindLUC_index = {'CO2': kLUC, 'BB': kLUC}

# other GHG
if (mod_kindGHG == 'one'):
    kind += ['non-CO2']
    kind_name += ['non-CO2']
    kind_color += ['#FFCC00']
    kCHI = kGHG + 1
elif (mod_kindGHG == 'RCP'):
    kind += ['CH4', 'N2O', 'HaloC']
    kind_name += ['Methane', 'Nitrous Oxide', 'Halocarbons']
    kind_color += ['#FF6600', '#FFCC00', '#FF9999']
    kCHI = kGHG + 3
else:
    kGHG = 0
    kCHI = max(kFF, kLUC) + 1

if (mod_kindGHG == 'RCP'):
    kindGHG_index = {'CH4': kGHG, 'N2O': kGHG + 1, 'HFC': kGHG + 2, 'PFC': kGHG + 2,
                     'ODS': kGHG + 2}
else:
    kindGHG_index = {'CH4': kGHG, 'N2O': kGHG, 'HFC': kGHG, 'PFC': kGHG, 'ODS': kGHG}

# active species
if (mod_kindCHI == 'one'):
    kind += ['OzPrec.']
    kind_name += ['Ozone Precursors']
    kind_color += ['#66FF66']
    kAER = kCHI + 1
elif (mod_kindCHI == 'all'):
    kind += ['NOx', 'CO', 'NMVOC']
    kind_name += ['Nitrogen Oxides', 'Carbon Monoxide',
                  'Non-Methane Volatile Organic Compounds']
    kind_color += ['#66FF66', '#009999', '#006600']
    kAER = kCHI + 3
else:
    kCHI = 0
    kAER = max(kFF, kLUC, kGHG) + 1

if (mod_kindCHI == 'all'):
    kindCHI_index = {'NOX': kCHI, 'CO': kCHI + 1, 'VOC': kCHI + 2}
else:
    kindCHI_index = {'NOX': kCHI, 'CO': kCHI, 'VOC': kCHI}

# aerosols
if (mod_kindAER == 'one'):
    kind += ['AER']
    kind_name += ['Aerosols']
    kind_color += ['#0000FF']
    kRF = kAER + 1
elif (mod_kindAER == 'all'):
    kind += ['SO2', 'NH3', 'OC', 'BC']
    kind_name += ['Sulfur Dioxide', 'Ammonia', 'Organic Carbon', 'Black Carbon']
    kind_color += ['#00CCFF', '#0000FF', '#660099', '#CC0066']
    kRF = kAER + 4
else:
    kAER = 0
    kRF = max(kFF, kLUC, kGHG, kCHI) + 1

if (mod_kindAER == 'all'):
    kindAER_index = {'SO2': kAER, 'NH3': kAER + 1, 'OC': kAER + 2, 'BC': kAER + 3}
else:
    kindAER_index = {'SO2': kAER, 'NH3': kAER, 'OC': kAER, 'BC': kAER}

# other radiative forcings
if (mod_kindRF == 'one'):
    kind += ['RFother']
    kind_name += ['Other RF']
    kind_color += ['#999999']
    kGE = kRF + 1
elif (mod_kindRF == 'two'):
    kind += ['RFant', 'RFnat']
    kind_name += ['Anthropogenic RF', 'Natural RF']
    kind_color += ['', '']
    kGE = kRF + 2
elif (mod_kindRF == 'all'):
    kind += ['RFcon', 'RFsol', 'RFvol']
    kind_name += ['Contrails RF', 'Solar RF', 'Volcanoes RF']
    kind_color += ['', '', '']
    kGE = kRF + 3
else:
    kRF = 0
    kGE = max(kFF, kLUC, kGHG, kCHI, kAER) + 1

if (mod_kindRF == 'all'):
    kindRF_index = {'RFcon': kRF, 'RFsol': kRF + 1, 'RFvol': kRF + 2}
elif (mod_kindRF == 'two'):
    kindRF_index = {'RFcon': kRF, 'RFsol': kRF + 1, 'RFvol': kRF + 1}
else:
    kindRF_index = {'RFcon': kRF, 'RFsol': kRF, 'RFvol': kRF}

# geoengineering
if (mod_kindGE == 'PUP'):
    kind += ['AFO', 'CCS', 'ALB', 'AER']
    kind_name += ['Aforestation', 'Carbon Capture and Storage', 'Surface Albedo',
                  'Sulfate Aerosols']
    kind_color += ['', '', '', '']
else:
    kGE = 0

nb_kind = len(kind)

# ===========
# A.4. Biomes
# ===========

if (mod_biomeSHR == 'w/FOR') & (mod_biomeURB == 'w/DES'):
    biome = ['DES+', 'FOR+', 'GRA', 'CRO', 'PAS']
    biome_name = ['Desert & Urban', 'Forest & Shrubland', 'Grassland', 'Cropland',
                  'Pasture']
    biome_color = ['', '', '', '', '']
    biome_index = {'des': 0, 'for': 1, 'shr': 1, 'gra': 2, 'cro': 3, 'pas': 4, 'urb': 0}
elif (mod_biomeSHR == 'w/FOR') & (mod_biomeURB == 'URB'):
    biome = ['DES', 'FOR+', 'GRA', 'CRO', 'PAS', 'URB']
    biome_name = ['Desert', 'Forest & Shrubland', 'Grassland', 'Cropland', 'Pasture',
                  'Urban']
    biome_color = ['', '', '', '', '', '']
    biome_index = {'des': 0, 'for': 1, 'shr': 1, 'gra': 2, 'cro': 3, 'pas': 4, 'urb': 5}
elif (mod_biomeSHR == 'w/GRA') & (mod_biomeURB == 'w/DES'):
    biome = ['DES+', 'FOR', 'GRA+', 'CRO', 'PAS']
    biome_name = ['Desert & Urban', 'Forest', 'Grassland & Shrubland', 'Cropland',
                  'Pasture']
    biome_color = ['', '', '', '', '']
    biome_index = {'des': 0, 'for': 1, 'shr': 2, 'gra': 2, 'cro': 3, 'pas': 4, 'urb': 0}
elif (mod_biomeSHR == 'w/GRA') & (mod_biomeURB == 'URB'):
    biome = ['DES', 'FOR', 'GRA+', 'CRO', 'PAS', 'URB']
    biome_name = ['Desert', 'Forest', 'Grassland & Shrubland', 'Cropland', 'Pasture',
                  'Urban']
    biome_color = ['', '', '', '', '', '']
    biome_index = {'des': 0, 'for': 1, 'shr': 2, 'gra': 2, 'cro': 3, 'pas': 4, 'urb': 5}
elif (mod_biomeSHR == 'SHR') & (mod_biomeURB == 'w/DES'):
    biome = ['DES+', 'FOR', 'SHR', 'GRA', 'CRO', 'PAS']
    biome_name = ['Desert & Urban', 'Forest', 'Shrubland', 'Grassland', 'Cropland',
                  'Pasture']
    biome_color = ['', '', '', '', '', '']
    biome_index = {'des': 0, 'for': 1, 'shr': 2, 'gra': 3, 'cro': 4, 'pas': 5, 'urb': 0}
elif (mod_biomeSHR == 'SHR') & (mod_biomeURB == 'URB'):
    biome = ['DES', 'FOR', 'SHR', 'GRA', 'CRO', 'PAS', 'URB']
    biome_name = ['Desert', 'Forest', 'Shrubland', 'Grassland', 'Cropland', 'Pasture',
                  'Urban']
    biome_color = ['', '', '', '', '', '', '']
    biome_index = {'des': 0, 'for': 1, 'shr': 2, 'gra': 3, 'cro': 4, 'pas': 5, 'urb': 6}
else:
    biome = ['all']
    biome_name = biome
    biome_color = []
    biome_index = {'des': 0, 'for': 0, 'shr': 0, 'gra': 0, 'cro': 0, 'pas': 0, 'urb': 0}

## for OSCAR v3 parameters
if mod_biomeV3:
    biome = ['FOR', 'GRA+', 'CRO', 'PAS', 'URB']
    biome_name = ['Forest', 'Non-Forest', 'Cropland', 'Pasture', 'Urban']
    biome_color = ['', '', '', '', '']
    biome_index = {'des': 1, 'for': 0, 'shr': 1, 'gra': 1, 'cro': 2, 'pas': 3, 'urb': 4}

nb_biome = len(biome)

# =========
# A.5. Halo
# =========

HFC = ['HFC23', 'HFC32', 'HFC125', 'HFC134a', 'HFC143a', 'HFC152a', 'HFC227ea',
       'HFC236fa', 'HFC245fa', 'HFC365mfc', 'HFC4310mee']
PFC = ['SF6', 'NF3', 'CF4', 'C2F6', 'C3F8', 'cC4F8', 'C4F10', 'C5F12', 'C6F14', 'C7F16']
ODS = ['CFC11', 'CFC12', 'CFC113', 'CFC114', 'CFC115', 'CCl4', 'CH3CCl3', 'HCFC22',
       'HCFC141b', 'HCFC142b', 'Halon1211', 'Halon1202', 'Halon1301', 'Halon2402',
       'CH3Br', 'CH3Cl']

nb_HFC = len(HFC)
nb_PFC = len(PFC)
nb_ODS = len(ODS)

##################################################
#   1. GREENHOUSE GASES
##################################################

EFF = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
               dtype=dty)  # {GtC/yr}
ECH4 = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                dtype=dty)  # {TgC/yr}
EN2O = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                dtype=dty)  # {TgN/yr}

# ==========
# 1.1. CDIAC
# ==========

# load CDIAC emissions
# from [Boden et al., 2012]
ind_cdiac = 310
kin = ['sol', 'liq', 'gas', 'cem', 'fla']

# total emissions
EFFcdiac = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                    dtype=dty)
for k in range(len(kin)):
    TMP = np.array([line for line in csv.reader(open(
        'data/EFossil_CDIAC/#DATA.EFossil_CDIAC.1751-2010_114reg0.EFF_' + kin[k] + '.csv',
        'r'))], dtype=dty)
    for i in range(114 + 1):
        EFFcdiac[51:ind_cdiac + 1, regionJ_index[i], 0, kFF, regionI_index[i]] += TMP[
                                                                                  :ind_cdiac - 51 + 1,
                                                                                  i]

# distribution among fuels
# TODO

# ========
# 1.2. EPA
# ========

# load EPA emissions
# see [EPA, 2012]
ind_epa = 310
sec_epa = ['agr_ferm', 'agr_manu', 'agr_othr', 'agr_rice', 'agr_soil', 'ene_burn',
           'ene_coal', 'ene_comb', 'ene_ngos', 'ene_othr', 'ind_acid']
sec_epa1 = ['agr_ferm', 'agr_manu', 'agr_othr', 'agr_rice', 'agr_soil', 'ene_coal',
            'ene_comb', 'ene_ngos', 'ene_othr', 'ind_acid']

# CH4
ECH4epa = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
for s in range(len(sec_epa1)):
    if os.path.isfile('data/EMethane_EPA/#DATA.EMethane_EPA.1990-' + str(
            1700 + ind_epa) + '_114reg0.ECH4_' + sec_epa1[s] + '.csv'):
        TMP = np.array([line for line in csv.reader(open(
            'data/EMethane_EPA/#DATA.EMethane_EPA.1990-' + str(
                1700 + ind_epa) + '_114reg0.ECH4_' + sec_epa1[s] + '.csv', 'r'))],
                       dtype=dty)
        for i in range(114 + 1):
            ECH4epa[290:ind_epa + 1, regionJ_index[i], 0, kindGHG_index['CH4'],
            regionI_index[i]] += TMP[:ind_epa - 290 + 1, i]

# N2O
EN2Oepa = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
for s in range(len(sec_epa1)):
    if os.path.isfile('data/ENitrousOx_EPA/#DATA.ENitrousOx_EPA.1990-' + str(
            1700 + ind_epa) + '_114reg0.EN2O_' + sec_epa1[s] + '.csv'):
        TMP = np.array([line for line in csv.reader(open(
            'data/ENitrousOx_EPA/#DATA.ENitrousOx_EPA.1990-' + str(
                1700 + ind_epa) + '_114reg0.EN2O_' + sec_epa1[s] + '.csv', 'r'))],
                       dtype=dty)
        for i in range(114 + 1):
            EN2Oepa[290:ind_epa + 1, regionJ_index[i], 0, kindGHG_index['N2O'],
            regionI_index[i]] += TMP[:ind_epa - 290 + 1, i]

# ==========
# 1.3. EDGAR
# ==========

# load emissions from EDGAR v4.2
# see [JRC, 2011]
ind_edgar = 308
sec_ehyde = ['oo', 'fc', 'fp', 'bc', 'in', 'al', 'an', 'aw', 'lf', 'sb', 'df']
sec_accmip = ['ooo', 'ene', 'ind', 'tra', 'shp', 'air', 'dom', 'slv', 'agr', 'awb', 'wst',
              'for', 'gra']

# FF
EFFedgar = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                    dtype=dty)
for s in range(1, len(sec_ehyde) - 2):
    if os.path.isfile('data/ECarbon_EDGAR/#DATA.ECarbon_EDGAR.1970-' + str(
            1700 + ind_edgar) + '_114reg0.ECO2_' + sec_ehyde[s] + '.csv'):
        TMP = np.array([line for line in csv.reader(open(
            'data/ECarbon_EDGAR/#DATA.ECarbon_EDGAR.1970-' + str(
                1700 + ind_edgar) + '_114reg0.ECO2_' + sec_ehyde[s] + '.csv', 'r'))],
                       dtype=dty)
        for i in range(114 + 1):
            EFFedgar[270:ind_edgar + 1, regionJ_index[i], 0, kFF,
            regionI_index[i]] += TMP[:ind_edgar - 270 + 1, i]

# CH4
ECH4edgar = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                     dtype=dty)
for s in range(1, len(sec_accmip) - 2):
    if os.path.isfile('data/EMethane_EDGAR/#DATA.EMethane_EDGAR.1970-' + str(
            1700 + ind_edgar) + '_114reg0.ECH4_' + sec_accmip[s] + '.csv'):
        TMP = np.array([line for line in csv.reader(open(
            'data/EMethane_EDGAR/#DATA.EMethane_EDGAR.1970-' + str(
                1700 + ind_edgar) + '_114reg0.ECH4_' + sec_accmip[s] + '.csv', 'r'))],
                       dtype=dty)
        for i in range(114 + 1):
            ECH4edgar[270:ind_edgar + 1, regionJ_index[i], 0, kindGHG_index['CH4'],
            regionI_index[i]] += TMP[:ind_edgar - 270 + 1, i]

# N2O
EN2Oedgar = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                     dtype=dty)
for s in range(1, len(sec_ehyde) - 2):
    if os.path.isfile('data/ENitrousOx_EDGAR/#DATA.ENitrousOx_EDGAR.1970-' + str(
            1700 + ind_edgar) + '_114reg0.EN2O_' + sec_ehyde[s] + '.csv'):
        TMP = np.array([line for line in csv.reader(open(
            'data/ENitrousOx_EDGAR/#DATA.ENitrousOx_EDGAR.1970-' + str(
                1700 + ind_edgar) + '_114reg0.EN2O_' + sec_ehyde[s] + '.csv', 'r'))],
                       dtype=dty)
        for i in range(114 + 1):
            EN2Oedgar[270:ind_edgar + 1, regionJ_index[i], 0, kindGHG_index['N2O'],
            regionI_index[i]] += TMP[:ind_edgar - 270 + 1, i]

# =============
# 1.4. EDGAR-FT
# =============

# load emissions from EDGAR-FT v4.2-FT2010
# see [JRC, 2013]

# FF
EFFeft = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
for s in range(1, len(sec_ehyde) - 2):
    if os.path.isfile(
            'data/ECarbon_EDGAR-FT/#DATA.ECarbon_EDGAR-FT.2008-2010_114reg0.ECO2_' +
            sec_ehyde[s] + '.csv'):
        TMP = np.array([line for line in csv.reader(open(
            'data/ECarbon_EDGAR-FT/#DATA.ECarbon_EDGAR-FT.2008-2010_114reg0.ECO2_' +
            sec_ehyde[s] + '.csv', 'r'))], dtype=dty)
        for i in range(114 + 1):
            EFFeft[308:310 + 1, regionJ_index[i], 0, kFF, regionI_index[i]] += TMP[
                                                                               :310 - 308 + 1,
                                                                               i]

# CH4
ECH4eft = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
for s in range(1, len(sec_accmip) - 2):
    if os.path.isfile(
            'data/EMethane_EDGAR-FT/#DATA.EMethane_EDGAR-FT.2008-2010_114reg0.ECH4_' +
            sec_accmip[s] + '.csv'):
        TMP = np.array([line for line in csv.reader(open(
            'data/EMethane_EDGAR-FT/#DATA.EMethane_EDGAR-FT.2008-2010_114reg0.ECH4_' +
            sec_accmip[s] + '.csv', 'r'))], dtype=dty)
        for i in range(114 + 1):
            ECH4eft[308:310 + 1, regionJ_index[i], 0, kindGHG_index['CH4'],
            regionI_index[i]] += TMP[:310 - 308 + 1, i]

# N2O
EN2Oeft = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
for s in range(1, len(sec_ehyde) - 2):
    if os.path.isfile(
            'data/ENitrousOx_EDGAR-FT/#DATA.ENitrousOx_EDGAR-FT.2008-2010_114reg0.EN2O_' +
            sec_ehyde[s] + '.csv'):
        TMP = np.array([line for line in csv.reader(open(
            'data/ENitrousOx_EDGAR-FT/#DATA.ENitrousOx_EDGAR-FT.2008-2010_114reg0.EN2O_' +
            sec_ehyde[s] + '.csv', 'r'))], dtype=dty)
        for i in range(114 + 1):
            EN2Oeft[308:310 + 1, regionJ_index[i], 0, kindGHG_index['N2O'],
            regionI_index[i]] += TMP[:310 - 308 + 1, i]

# ===========
# 1.5. ACCMIP
# ===========

# load emissions from ACCMIP
# from [Lamarque et al., 2010]

# CH4
ECH4accmip = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                      dtype=dty)
p_ECH4_bio = np.zeros([nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
for s in range(1, len(sec_accmip) - 2):
    if os.path.isfile(
            'data/EMethane_ACCMIP/#DATA.EMethane_ACCMIP.1850-2000_114reg0.ECH4_' +
            sec_accmip[s] + '.csv'):
        TMP = np.array([line for line in csv.reader(open(
            'data/EMethane_ACCMIP/#DATA.EMethane_ACCMIP.1850-2000_114reg0.ECH4_' +
            sec_accmip[s] + '.csv', 'r'))], dtype=dty)
        for i in range(114 + 1):
            ECH4accmip[150:300 + 1, regionJ_index[i], 0, kindGHG_index['CH4'],
            regionI_index[i]] += TMP[:300 - 150 + 1, i]
            if sec_accmip[s] in ['agr', 'awb', 'wst']:
                p_ECH4_bio[regionJ_index[i], 0, kindGHG_index['CH4'], regionI_index[i]] += \
                TMP[0, i]
p_ECH4_bio /= ECH4accmip[150]
p_ECH4_bio[np.isnan(p_ECH4_bio) | np.isinf(p_ECH4_bio)] = 0

# ===============
# 1.6. EDGAR-HYDE
# ===============

# load emissions from EDGAR-HYDE v1.3
# from [van Aardenne et al., 2001]

# N2O
EN2Oehyde = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                     dtype=dty)
p_EN2O_bio = np.zeros([nb_regionJ, nb_sector, nb_kind, nb_regionI], dtype=dty)
for s in range(1, len(sec_ehyde) - 2):
    if os.path.isfile(
            'data/ENitrousOx_EDGAR-HYDE/#DATA.ENitrousOx_EDGAR-HYDE.1890-1990_114reg0.EN2O_' +
            sec_ehyde[s] + '.csv'):
        TMP = np.array([line for line in csv.reader(open(
            'data/ENitrousOx_EDGAR-HYDE/#DATA.ENitrousOx_EDGAR-HYDE.1890-1990_114reg0.EN2O_' +
            sec_ehyde[s] + '.csv', 'r'))], dtype=dty)
        for i in range(114 + 1):
            EN2Oehyde[190:290 + 1, regionJ_index[i], 0, kindGHG_index['N2O'],
            regionI_index[i]] += TMP[:290 - 190 + 1, i]
            if sec_ehyde[s] in ['al', 'an', 'aw', 'lf']:
                p_EN2O_bio[regionJ_index[i], 0, kindGHG_index['N2O'], regionI_index[i]] += \
                TMP[0, i]
p_EN2O_bio /= EN2Oehyde[190]
p_EN2O_bio[np.isnan(p_EN2O_bio) | np.isinf(p_EN2O_bio)] = 0

# ==============
# 1.7. Stern1998
# ==============

# load emissions from [Stern et al., 1998]
ECH4stern = np.zeros([ind_cdiac + 1], dtype=dty)
TMP = np.array([line for line in csv.reader(
    open('data/EMethane_Stern1998/#DATA.EMethane_Stern1998.1860-1994_(7sec).ECH4.csv',
         'r'))][1:], dtype=dty)
lgd = [line for line in csv.reader(
    open('data/EMethane_Stern1998/#DATA.EMethane_Stern1998.1860-1994_(7sec).ECH4.csv',
         'r'))][0]
for s in range(len(lgd)):
    if not lgd[s] in ['Biomass Burning']:
        ECH4stern[160:294 + 1] += TMP[:294 - 160 + 1, s]

# =================
# 1.8. Davidson2009
# =================

# load emissions from [Davidson et al., 2009]
EN2Odavidson = np.zeros([ind_cdiac + 1], dtype=dty)
TMP = np.array([line for line in csv.reader(open(
    'data/ENitrousOx_Davidson2009/#DATA.ENitrousOx_Davidson2009.1860-2005_(5sec).EN2O.csv',
    'r'))][1:], dtype=dty)
lgd = [line for line in csv.reader(open(
    'data/ENitrousOx_Davidson2009/#DATA.ENitrousOx_Davidson2009.1860-2005_(5sec).EN2O.csv',
    'r'))][0]
for s in range(len(lgd)):
    if lgd[s] in ['nyl_prod', 'ff_burn', 'biogen']:
        EN2Odavidson[160:305 + 1] += TMP[:305 - 160 + 1, s]

# global rescaling of EDGAR-HYDE emissions
EN2OehydeR = EN2Oehyde.copy()
if True:
    EN2OehydeR[190:290 + 1, ...] *= (EN2Odavidson[190:290 + 1] / np.sum(
        np.sum(np.sum(np.sum(EN2Oehyde[190:290 + 1, ...], 4), 3), 2), 1))[:, np.newaxis,
                                    np.newaxis, np.newaxis, np.newaxis]

# =========
# 1.9. SRES
# =========

# initialization of projected drivers
for VAR in ['FF', 'CH4', 'N2O']:
    exec(
        'E' + VAR + 'proj = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)')

# projection of emissions under SRES scenarios
# from [IPCC, 2000]

# FF
if (scen_EFF[:4] == 'SRES') & (ind_final > ind_cdiac):
    TMP = np.array([line for line in csv.reader(open(
        'data/EFossil_SRES/#DATA.EFossil_SRES.2000-2100_4reg0.' + scen_EFF[
                                                                  5:] + '_EFF.csv',
        'r'))], dtype=dty)
    for i in range(4 + 1):
        if (mod_regionI == 'SRES4') & (mod_regionJ == 'SRES4'):
            EFFproj[300:min(ind_final, 400) + 1, i, 0, kFF, i] += TMP[:min(ind_final,
                                                                           400) - 300 + 1,
                                                                  i]
        elif (mod_regionI == 'SRES4') & (mod_regionJ != 'SRES4'):
            EFFproj[300:min(ind_final, 400) + 1, 0, 0, kFF, i] += TMP[:min(ind_final,
                                                                           400) - 300 + 1,
                                                                  i]
        elif (mod_regionI != 'SRES4') & (mod_regionJ == 'SRES4'):
            EFFproj[300:min(ind_final, 400) + 1, i, 0, kFF, 0] += TMP[:min(ind_final,
                                                                           400) - 300 + 1,
                                                                  i]
        elif (mod_regionI != 'SRES4') & (mod_regionJ != 'SRES4'):
            EFFproj[300:min(ind_final, 400) + 1, 0, 0, kFF, 0] += TMP[:min(ind_final,
                                                                           400) - 300 + 1,
                                                                  i]

# CH4
if (scen_ECH4[:4] == 'SRES') & (ind_final > ind_cdiac):
    TMP = np.array([line for line in csv.reader(open(
        'data/EMethane_SRES/#DATA.EMethane_SRES.2000-2100_4reg0.' + scen_ECH4[
                                                                    5:] + '_ECH4.csv',
        'r'))], dtype=dty)
    for i in range(4 + 1):
        if (mod_regionI == 'SRES4') & (mod_regionJ == 'SRES4'):
            ECH4proj[300:min(ind_final, 400) + 1, i, 0, kindGHG_index['CH4'], i] += TMP[
                                                                                    :min(
                                                                                        ind_final,
                                                                                        400) - 300 + 1,
                                                                                    i]
        elif (mod_regionI == 'SRES4') & (mod_regionJ != 'SRES4'):
            ECH4proj[300:min(ind_final, 400) + 1, 0, 0, kindGHG_index['CH4'], i] += TMP[
                                                                                    :min(
                                                                                        ind_final,
                                                                                        400) - 300 + 1,
                                                                                    i]
        elif (mod_regionI != 'SRES4') & (mod_regionJ == 'SRES4'):
            ECH4proj[300:min(ind_final, 400) + 1, i, 0, kindGHG_index['CH4'], 0] += TMP[
                                                                                    :min(
                                                                                        ind_final,
                                                                                        400) - 300 + 1,
                                                                                    i]
        elif (mod_regionI != 'SRES4') & (mod_regionJ != 'SRES4'):
            ECH4proj[300:min(ind_final, 400) + 1, 0, 0, kindGHG_index['CH4'], 0] += TMP[
                                                                                    :min(
                                                                                        ind_final,
                                                                                        400) - 300 + 1,
                                                                                    i]

# N2O
if (scen_EN2O[:4] == 'SRES') & (ind_final > ind_cdiac):
    TMP = np.array([line for line in csv.reader(open(
        'data/ENitrousOx_SRES/#DATA.ENitrousOx_SRES.2000-2100_4reg0.' + scen_EN2O[
                                                                        5:] + '_EN2O.csv',
        'r'))], dtype=dty)
    for i in range(4 + 1):
        if (mod_regionI == 'SRES4') & (mod_regionJ == 'SRES4'):
            EN2Oproj[300:min(ind_final, 400) + 1, i, 0, kindGHG_index['N2O'], i] += TMP[
                                                                                    :min(
                                                                                        ind_final,
                                                                                        400) - 300 + 1,
                                                                                    i]
        elif (mod_regionI == 'SRES4') & (mod_regionJ != 'SRES4'):
            EN2Oproj[300:min(ind_final, 400) + 1, 0, 0, kindGHG_index['N2O'], i] += TMP[
                                                                                    :min(
                                                                                        ind_final,
                                                                                        400) - 300 + 1,
                                                                                    i]
        elif (mod_regionI != 'SRES4') & (mod_regionJ == 'SRES4'):
            EN2Oproj[300:min(ind_final, 400) + 1, i, 0, kindGHG_index['N2O'], 0] += TMP[
                                                                                    :min(
                                                                                        ind_final,
                                                                                        400) - 300 + 1,
                                                                                    i]
        elif (mod_regionI != 'SRES4') & (mod_regionJ != 'SRES4'):
            EN2Oproj[300:min(ind_final, 400) + 1, 0, 0, kindGHG_index['N2O'], 0] += TMP[
                                                                                    :min(
                                                                                        ind_final,
                                                                                        400) - 300 + 1,
                                                                                    i]

# =========
# 1.10. RCP
# =========

# projection of emissions under RCP scenarios
# from [Meinshausen et al., 2011]

# FF
if (scen_EFF[:3] == 'RCP') & (ind_final > ind_cdiac):
    TMP = np.array([line for line in csv.reader(open(
        'data/EFossil_RCP/#DATA.EFossil_RCP.2000-2100_5reg0.rcp' + scen_EFF[3] + scen_EFF[
            5] + '_EFF.csv', 'r'))], dtype=dty)
    for i in range(5 + 1):
        if (mod_regionI == 'RCP5') & (mod_regionJ == 'RCP5'):
            EFFproj[300:min(ind_final, 400) + 1, i, 0, kFF, i] += TMP[:min(ind_final,
                                                                           400) - 300 + 1,
                                                                  i]
        elif (mod_regionI == 'RCP5') & (mod_regionJ != 'RCP5'):
            EFFproj[300:min(ind_final, 400) + 1, 0, 0, kFF, i] += TMP[:min(ind_final,
                                                                           400) - 300 + 1,
                                                                  i]
        elif (mod_regionI != 'RCP5') & (mod_regionJ == 'RCP5'):
            EFFproj[300:min(ind_final, 400) + 1, i, 0, kFF, 0] += TMP[:min(ind_final,
                                                                           400) - 300 + 1,
                                                                  i]
        elif (mod_regionI != 'RCP5') & (mod_regionJ != 'RCP5'):
            EFFproj[300:min(ind_final, 400) + 1, 0, 0, kFF, 0] += TMP[:min(ind_final,
                                                                           400) - 300 + 1,
                                                                  i]

# CH4
if (scen_ECH4[:3] == 'RCP') & (ind_final > ind_cdiac):
    for s in range(1, len(sec_accmip) - 2):
        if os.path.isfile(
                'data/EMethane_RCP/#DATA.EMethane_RCP.2000-2100_5reg0.rcp' + scen_ECH4[
                    3] + scen_ECH4[5] + '_ECH4_' + sec_accmip[s] + '.csv'):
            TMP = np.array([line for line in csv.reader(open(
                'data/EMethane_RCP/#DATA.EMethane_RCP.2000-2100_5reg0.rcp' + scen_ECH4[
                    3] + scen_ECH4[5] + '_ECH4_' + sec_accmip[s] + '.csv', 'r'))],
                           dtype=dty)
            for i in range(5 + 1):
                if (mod_regionI == 'RCP5') & (mod_regionJ == 'RCP5'):
                    ECH4proj[300:min(ind_final, 400) + 1, i, 0, kindGHG_index['CH4'],
                    i] += TMP[:min(ind_final, 400) - 300 + 1, i]
                elif (mod_regionI == 'RCP5') & (mod_regionJ != 'RCP5'):
                    ECH4proj[300:min(ind_final, 400) + 1, 0, 0, kindGHG_index['CH4'],
                    i] += TMP[:min(ind_final, 400) - 300 + 1, i]
                elif (mod_regionI != 'RCP5') & (mod_regionJ == 'RCP5'):
                    ECH4proj[300:min(ind_final, 400) + 1, i, 0, kindGHG_index['CH4'],
                    0] += TMP[:min(ind_final, 400) - 300 + 1, i]
                elif (mod_regionI != 'RCP5') & (mod_regionJ != 'RCP5'):
                    ECH4proj[300:min(ind_final, 400) + 1, 0, 0, kindGHG_index['CH4'],
                    0] += TMP[:min(ind_final, 400) - 300 + 1, i]

# N2O
if (scen_EN2O[:3] == 'RCP') & (ind_final > ind_cdiac):
    TMP = np.array([line for line in csv.reader(open(
        'data/ENitrousOx_RCP/#DATA.ENitrousOx_RCP.2000-2100_5reg0.rcp' + scen_EN2O[3] +
        scen_EN2O[5] + '_EN2O.csv', 'r'))], dtype=dty)
    for i in range(5 + 1):
        if (mod_regionI == 'RCP5') & (mod_regionJ == 'RCP5'):
            EN2Oproj[300:min(ind_final, 400) + 1, i, 0, kindGHG_index['N2O'], i] += TMP[
                                                                                    :min(
                                                                                        ind_final,
                                                                                        400) - 300 + 1,
                                                                                    i]
        elif (mod_regionI == 'RCP5') & (mod_regionJ != 'RCP5'):
            EN2Oproj[300:min(ind_final, 400) + 1, 0, 0, kindGHG_index['N2O'], i] += TMP[
                                                                                    :min(
                                                                                        ind_final,
                                                                                        400) - 300 + 1,
                                                                                    i]
        elif (mod_regionI != 'RCP5') & (mod_regionJ == 'RCP5'):
            EN2Oproj[300:min(ind_final, 400) + 1, i, 0, kindGHG_index['N2O'], 0] += TMP[
                                                                                    :min(
                                                                                        ind_final,
                                                                                        400) - 300 + 1,
                                                                                    i]
        elif (mod_regionI != 'RCP5') & (mod_regionJ != 'RCP5'):
            EN2Oproj[300:min(ind_final, 400) + 1, 0, 0, kindGHG_index['N2O'], 0] += TMP[
                                                                                    :min(
                                                                                        ind_final,
                                                                                        400) - 300 + 1,
                                                                                    i]

# =================
# 1.A. PAST DATASET
# =================

# datasets mixed following trends
for VAR in ['FF', 'CH4', 'N2O']:

    if VAR in ['FF']:

        # with CDIAC as reference
        if (data_EFF == 'CDIAC'):
            exec('E' + VAR + 'past = E' + VAR + 'cdiac.copy()')

        # with EDGAR as reference
        elif (data_EFF == 'EDGAR'):
            exec('E' + VAR + 'past = E' + VAR + 'edgar.copy()')
            # follow EDGAR-FT variations after 2008
            for t in range(ind_edgar + 1, ind_cdiac + 1):
                exec(
                    'E' + VAR + 'past[t,...] = E' + VAR + 'past[t-1,...] * E' + VAR + 'eft[t,...]/E' + VAR + 'eft[t-1,...]')
                exec(
                    'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
                exec(
                    'E' + VAR + 'past[t,...] *= np.sum(E' + VAR + 'past[t-1,...])/np.sum(E' + VAR + 'past[t,...]) * np.sum(E' + VAR + 'eft[t,...])/np.sum(E' + VAR + 'eft[t-1,...])')
            exec(
                'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
            # follow CDIAC variations before 1970
            for t in range(50, 270)[::-1]:
                exec(
                    'E' + VAR + 'past[t,...] = E' + VAR + 'past[t+1,...] * E' + VAR + 'cdiac[t,...]/E' + VAR + 'cdiac[t+1,...]')
                exec(
                    'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
                exec(
                    'E' + VAR + 'past[t,...] *= np.sum(E' + VAR + 'past[t+1,...])/np.sum(E' + VAR + 'past[t,...]) * np.sum(E' + VAR + 'cdiac[t,...])/np.sum(E' + VAR + 'cdiac[t+1,...])')
            exec(
                'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')

    elif VAR in ['CH4']:

        # with EDGAR as reference
        if (data_ECH4 == 'EDGAR'):
            exec('E' + VAR + 'past = E' + VAR + 'edgar.copy()')
            # follow EDGAR-FT variations after 2008
            for t in range(ind_edgar + 1, ind_cdiac + 1):
                exec(
                    'E' + VAR + 'past[t,...] = E' + VAR + 'past[t-1,...] * E' + VAR + 'eft[t,...]/E' + VAR + 'eft[t-1,...]')
                exec(
                    'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
                exec(
                    'E' + VAR + 'past[t,...] *= np.sum(E' + VAR + 'past[t-1,...])/np.sum(E' + VAR + 'past[t,...]) * np.sum(E' + VAR + 'eft[t,...])/np.sum(E' + VAR + 'eft[t-1,...])')
            exec(
                'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
            # follow ACCMIP variations before 1970
            for t in range(150, 270)[::-1]:
                exec(
                    'E' + VAR + 'past[t,...] = E' + VAR + 'past[t+1,...] * E' + VAR + 'accmip[t,...]/E' + VAR + 'accmip[t+1,...]')
                exec(
                    'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
                exec(
                    'E' + VAR + 'past[t,...] *= np.sum(E' + VAR + 'past[t+1,...])/np.sum(E' + VAR + 'past[t,...]) * np.sum(E' + VAR + 'accmip[t,...])/np.sum(E' + VAR + 'accmip[t+1,...])')
            exec(
                'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
            # linear extrapolation before 1850
            for t in range(50, 150):
                exec(
                    'E' + VAR + 'past[t,...] = (1-p_E' + VAR + '_bio) * E' + VAR + 'past[150,...] * (t-50)/float(150-50)')
            for t in range(0, 150):
                exec(
                    'E' + VAR + 'past[t,...] += p_E' + VAR + '_bio * E' + VAR + 'past[150,...] * (t--200)/float(150--200)')
            exec('E' + VAR + '_0 = E' + VAR + 'past[50,...]')

        # with ACCMIP as reference
        elif (data_ECH4 == 'ACCMIP'):
            exec('E' + VAR + 'past = E' + VAR + 'accmip.copy()')
            # follow EDGAR variations after 2000
            for t in range(300 + 1, ind_edgar + 1):
                exec(
                    'E' + VAR + 'past[t,...] = E' + VAR + 'past[t-1,...] * E' + VAR + 'edgar[t,...]/E' + VAR + 'edgar[t-1,...]')
                exec(
                    'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
                exec(
                    'E' + VAR + 'past[t,...] *= np.sum(E' + VAR + 'past[t-1,...])/np.sum(E' + VAR + 'past[t,...]) * np.sum(E' + VAR + 'edgar[t,...])/np.sum(E' + VAR + 'edgar[t-1,...])')
            exec(
                'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
            # follow EDGAR-FT variations after 2008
            for t in range(ind_edgar + 1, ind_cdiac + 1):
                exec(
                    'E' + VAR + 'past[t,...] = E' + VAR + 'past[t-1,...] * E' + VAR + 'eft[t,...]/E' + VAR + 'eft[t-1,...]')
                exec(
                    'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
                exec(
                    'E' + VAR + 'past[t,...] *= np.sum(E' + VAR + 'past[t-1,...])/np.sum(E' + VAR + 'past[t,...]) * np.sum(E' + VAR + 'eft[t,...])/np.sum(E' + VAR + 'eft[t-1,...])')
            exec(
                'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
            # linear extrapolation before 1850
            for t in range(50, 150):
                exec(
                    'E' + VAR + 'past[t,...] = (1-p_E' + VAR + '_bio) * E' + VAR + 'past[150,...] * (t-50)/float(150-50)')
            for t in range(0, 150):
                exec(
                    'E' + VAR + 'past[t,...] += p_E' + VAR + '_bio * E' + VAR + 'past[150,...] * (t--200)/float(150--200)')
            exec('E' + VAR + '_0 = E' + VAR + 'past[50,...]')

        # with EPA as reference
        elif (data_ECH4 == 'EPA'):
            exec('E' + VAR + 'past = E' + VAR + 'epa.copy()')
            # follow ACCMIP variations before 1990
            for t in range(150, 290)[::-1]:
                exec(
                    'E' + VAR + 'past[t,...] = E' + VAR + 'past[t+1,...] * E' + VAR + 'accmip[t,...]/E' + VAR + 'accmip[t+1,...]')
                exec(
                    'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
                exec(
                    'E' + VAR + 'past[t,...] *= np.sum(E' + VAR + 'past[t+1,...])/np.sum(E' + VAR + 'past[t,...]) * np.sum(E' + VAR + 'accmip[t,...])/np.sum(E' + VAR + 'accmip[t+1,...])')
            exec(
                'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
            # linear extrapolation before 1850
            for t in range(50, 150):
                exec(
                    'E' + VAR + 'past[t,...] = (1-p_E' + VAR + '_bio) * E' + VAR + 'past[150,...] * (t-50)/float(150-50)')
            for t in range(0, 150):
                exec(
                    'E' + VAR + 'past[t,...] += p_E' + VAR + '_bio * E' + VAR + 'past[150,...] * (t--200)/float(150--200)')
            exec('E' + VAR + '_0 = E' + VAR + 'past[50,...]')

    elif VAR in ['N2O']:

        # with EDGAR as reference
        if (data_EN2O == 'EDGAR'):
            exec('E' + VAR + 'past = E' + VAR + 'edgar.copy()')
            # follow EDGAR-FT variations after 2008
            for t in range(ind_edgar + 1, ind_cdiac + 1):
                exec(
                    'E' + VAR + 'past[t,...] = E' + VAR + 'past[t-1,...] * E' + VAR + 'eft[t,...]/E' + VAR + 'eft[t-1,...]')
                exec(
                    'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
                exec(
                    'E' + VAR + 'past[t,...] *= np.sum(E' + VAR + 'past[t-1,...])/np.sum(E' + VAR + 'past[t,...]) * np.sum(E' + VAR + 'eft[t,...])/np.sum(E' + VAR + 'eft[t-1,...])')
            exec(
                'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
            # follow EDGAR-HYDE variations before 1970
            for t in range(190, 270)[::-1]:
                exec(
                    'E' + VAR + 'past[t,...] = E' + VAR + 'past[t+1,...] * E' + VAR + 'ehydeR[t,...]/E' + VAR + 'ehydeR[t+1,...]')
                exec(
                    'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
                exec(
                    'E' + VAR + 'past[t,...] *= np.sum(E' + VAR + 'past[t+1,...])/np.sum(E' + VAR + 'past[t,...]) * np.sum(E' + VAR + 'ehydeR[t,...])/np.sum(E' + VAR + 'ehydeR[t+1,...])')
            exec(
                'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
            # linear extrapolation before 1890
            for t in range(50, 190):
                exec(
                    'E' + VAR + 'past[t,...] = (1-p_E' + VAR + '_bio) * E' + VAR + 'past[190,...] * (t-50)/float(190-50)')
            for t in range(0, 190):
                exec(
                    'E' + VAR + 'past[t,...] += p_E' + VAR + '_bio * E' + VAR + 'past[190,...] * (t--200)/float(190--200)')
            exec('E' + VAR + '_0 = E' + VAR + 'past[50,...]')

        # with EPA as reference
        elif (data_EN2O == 'EPA'):
            exec('E' + VAR + 'past = E' + VAR + 'epa.copy()')
            # follow EDGAR-HYDE variations before 1990
            for t in range(190, 290)[::-1]:
                exec(
                    'E' + VAR + 'past[t,...] = E' + VAR + 'past[t+1,...] * E' + VAR + 'ehydeR[t,...]/E' + VAR + 'ehydeR[t+1,...]')
                exec(
                    'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
                exec(
                    'E' + VAR + 'past[t,...] *= np.sum(E' + VAR + 'past[t+1,...])/np.sum(E' + VAR + 'past[t,...]) * np.sum(E' + VAR + 'ehydeR[t,...])/np.sum(E' + VAR + 'ehydeR[t+1,...])')
            exec(
                'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
            # linear extrapolation before 1850
            for t in range(50, 190):
                exec(
                    'E' + VAR + 'past[t,...] = (1-p_E' + VAR + '_bio) * E' + VAR + 'past[190,...] * (t-50)/float(190-50)')
            for t in range(0, 190):
                exec(
                    'E' + VAR + 'past[t,...] += p_E' + VAR + '_bio * E' + VAR + 'past[190,...] * (t--200)/float(190--200)')
            exec('E' + VAR + '_0 = E' + VAR + 'past[50,...]')

    # cut past dataset to right length
    exec(
        'E' + VAR + '[:min(ind_cdiac,ind_final)+1,...] = E' + VAR + 'past[:min(ind_cdiac,ind_final)+1,...]')

# ==================
# 1.B. FINAL DATASET
# ==================

# datasets mixed following various criteria
for VAR in ['FF', 'CH4', 'N2O']:
    exec('scen = scen_E' + VAR)

    # stop emissions
    if (scen == 'stop') & (ind_final > ind_cdiac):
        if VAR in ['FF']:
            exec('E' + VAR + '[ind_cdiac+1:,...] = 0')
        else:
            exec('E' + VAR + '[ind_cdiac+1:,...] = E' + VAR + '_0[np.newaxis,...]')

    # constant emissions
    elif (scen == 'cst') & (ind_final > ind_cdiac):
        exec(
            'E' + VAR + '[ind_cdiac+1:,...] = E' + VAR + '[ind_cdiac,...][np.newaxis,...]')

        # RCP or SRES scenarios
    elif ((scen[:4] == 'SRES') | (scen[:3] == 'RCP')) & (ind_final > ind_cdiac):

        # raw discontinuity
        if (mod_DATAscen == 'raw'):
            exec('E' + VAR + '[ind_cdiac+1:,...] = E' + VAR + 'proj[ind_cdiac+1:,...]')

        # offset at transition point
        elif (mod_DATAscen == 'offset'):
            exec(
                'E' + VAR + '[ind_cdiac+1:,...] = E' + VAR + 'proj[ind_cdiac+1:,...] - E' + VAR + 'proj[ind_cdiac,...] + E' + VAR + '[ind_cdiac,...]')
            for t in range(ind_cdiac + 1, ind_final + 1):
                exec('def_regI = bool(np.sum(E' + VAR + 'proj[t,:,...,1:]))')
                exec('def_regJ = bool(np.sum(E' + VAR + 'proj[t,1:,...,:]))')
                if not def_regI:
                    exec('E' + VAR + '[t,:,...,0] += np.sum(E' + VAR + '[t,:,...,1:],-1)')
                    exec('E' + VAR + '[t,:,...,1:] = 0')
                if not def_regJ:
                    exec('E' + VAR + '[t,0,...,:] += np.sum(E' + VAR + '[t,1:,...,:],0)')
                    exec('E' + VAR + '[t,1:,...,:] = 0')

                    # linear transition over N years
        elif (mod_DATAscen[:6] == 'smooth'):
            N = int(mod_DATAscen[6:])
            if (ind_final >= ind_cdiac + N):
                for t in range(ind_cdiac + 1, ind_cdiac + N):
                    exec(
                        'E' + VAR + '[t,...] = (1-(t-ind_cdiac)/float(N)) * E' + VAR + '[ind_cdiac,...] + (t-ind_cdiac)/float(N) * E' + VAR + 'proj[ind_cdiac+N,...]')
                    exec('def_regI = bool(np.sum(E' + VAR + 'proj[t,:,...,1:]))')
                    exec('def_regJ = bool(np.sum(E' + VAR + 'proj[t,1:,...,:]))')
                    if not def_regI:
                        exec(
                            'E' + VAR + '[t,:,...,0] += np.sum(E' + VAR + '[t,:,...,1:],-1)')
                        exec('E' + VAR + '[t,:,...,1:] = 0')
                    if not def_regJ:
                        exec(
                            'E' + VAR + '[t,0,...,:] += np.sum(E' + VAR + '[t,1:,...,:],0)')
                        exec('E' + VAR + '[t,1:,...,:] = 0')
                exec(
                    'E' + VAR + '[ind_cdiac+N:,...] = E' + VAR + 'proj[ind_cdiac+N:,...]')

        # follow trends of projection
        elif (mod_DATAscen == 'trends'):
            for t in range(ind_cdiac + 1, ind_final + 1):
                exec('def_regI = bool(np.sum(E' + VAR + 'proj[t,:,...,1:]))')
                exec('def_regJ = bool(np.sum(E' + VAR + 'proj[t,1:,...,:]))')
                if def_regI and def_regJ:
                    exec(
                        'E' + VAR + '[t,...] = E' + VAR + '[t-1,...] * E' + VAR + 'proj[t,...]/E' + VAR + 'proj[t-1,...]')
                    exec(
                        'E' + VAR + '[np.isnan(E' + VAR + ')|np.isinf(E' + VAR + ')] = 0')
                    exec(
                        'E' + VAR + '[t,...] *= np.sum(E' + VAR + '[t-1,...])/np.sum(E' + VAR + '[t,...]) * np.sum(E' + VAR + 'proj[t,...])/np.sum(E' + VAR + 'proj[t-1,...])')
                elif not def_regI and def_regJ:
                    exec(
                        'E' + VAR + '[t,:,...,0] = np.sum(E' + VAR + '[t-1,:,...,:],-1) * np.sum(E' + VAR + 'proj[t,:,...,:],-1)/np.sum(E' + VAR + 'proj[t-1,:,...,:],-1)')
                    exec(
                        'E' + VAR + '[np.isnan(E' + VAR + ')|np.isinf(E' + VAR + ')] = 0')
                    exec(
                        'E' + VAR + '[t,...] *= np.sum(E' + VAR + '[t-1,...])/np.sum(E' + VAR + '[t,...]) * np.sum(E' + VAR + 'proj[t,...])/np.sum(E' + VAR + 'proj[t-1,...])')
                elif def_regI and not def_regJ:
                    exec(
                        'E' + VAR + '[t,0,...,:] = np.sum(E' + VAR + '[t-1,:,...,:],0) * np.sum(E' + VAR + 'proj[t,:,...,:],0)/np.sum(E' + VAR + 'proj[t-1,:,...,:],0)')
                    exec(
                        'E' + VAR + '[np.isnan(E' + VAR + ')|np.isinf(E' + VAR + ')] = 0')
                    exec(
                        'E' + VAR + '[t,...] *= np.sum(E' + VAR + '[t-1,...])/np.sum(E' + VAR + '[t,...]) * np.sum(E' + VAR + 'proj[t,...])/np.sum(E' + VAR + 'proj[t-1,...])')
                elif not def_regI and not def_regJ:
                    exec(
                        'E' + VAR + '[t,0,...,0] = np.sum(np.sum(E' + VAR + '[t-1,:,...,:],-1),0) * np.sum(np.sum(E' + VAR + 'proj[t,:,...,:],-1),0)/np.sum(np.sum(E' + VAR + 'proj[t-1,:,...,:],-1),0)')
            exec('E' + VAR + '[np.isnan(E' + VAR + ')|np.isinf(E' + VAR + ')] = 0')

# delete individual datasets
for VAR in ['FF']:
    exec(
        'del E' + VAR + 'cdiac,E' + VAR + 'edgar,E' + VAR + 'eft,E' + VAR + 'past,E' + VAR + 'proj')
for VAR in ['CH4']:
    exec(
        'del E' + VAR + 'epa,E' + VAR + 'edgar,E' + VAR + 'eft,E' + VAR + 'accmip,E' + VAR + 'past,E' + VAR + 'proj,E' + VAR + 'stern,p_E' + VAR + '_bio')
for VAR in ['N2O']:
    exec(
        'del E' + VAR + 'epa,E' + VAR + 'edgar,E' + VAR + 'eft,E' + VAR + 'ehyde,E' + VAR + 'past,E' + VAR + 'proj,E' + VAR + 'davidson,p_E' + VAR + '_bio')

# ===========
# 1.11. PETERS
# ===========

# TODO


##################################################
#   2. LAND-USE CHANGE
##################################################

LUC = np.zeros(
    [ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_biome, nb_biome],
    dtype=dty)  # {Mha/yr}
HARV = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_biome],
                dtype=dty)  # {GtC/yr}
SHIFT = np.zeros(
    [ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_biome, nb_biome],
    dtype=dty)  # {Mha/yr}

# =========
# 2.1. LUH1
# =========

# load land-use data from LUH1
# from [Hurtt et al., 2011] and updated
bio = ['des', 'for', 'shr', 'gra', 'cro', 'pas', 'urb']

# LUC
if (data_LULCC[:3] == 'LUH'):
    for b1 in range(len(bio)):
        for b2 in range(len(bio)):
            if os.path.isfile(
                    'data/LandUse_' + data_LULCC + '/#DATA.LandUse_' + data_LULCC + '_' + mod_LSNKcover + '.1501-2015_114reg1.LUC_' +
                    bio[b1] + '2' + bio[b2] + '.csv'):
                TMP = np.array([line for line in csv.reader(open(
                    'data/LandUse_LUH1/#DATA.LandUse_' + data_LULCC + '_' + mod_LSNKcover + '.1501-2015_114reg1.LUC_' +
                    bio[b1] + '2' + bio[b2] + '.csv', 'r'))], dtype=dty)
                for i in range(1, 114 + 1):
                    LUC[1:min(ind_final, ind_cdiac) + 1, regionJ_index[i], 0, kLUC,
                    regionI_index[i], biome_index[bio[b1]], biome_index[bio[b2]]] += TMP[
                                                                                     200:200 + min(
                                                                                         ind_final,
                                                                                         ind_cdiac) - 1 + 1,
                                                                                     i - 1]

# HARV
if (data_LULCC[:3] == 'LUH'):
    for b in range(len(bio)):
        if os.path.isfile(
                'data/LandUse_' + data_LULCC + '/#DATA.LandUse_' + data_LULCC + '_' + mod_LSNKcover + '.1501-2015_114reg1.HARV_' +
                bio[b] + '.csv'):
            TMP = np.array([line for line in csv.reader(open(
                'data/LandUse_LUH1/#DATA.LandUse_' + data_LULCC + '_' + mod_LSNKcover + '.1501-2015_114reg1.HARV_' +
                bio[b] + '.csv', 'r'))], dtype=dty)
            for i in range(1, 114 + 1):
                HARV[1:min(ind_final, ind_cdiac) + 1, regionJ_index[i], 0, kLUC,
                regionI_index[i], biome_index[bio[b]]] += TMP[200:200 + min(ind_final,
                                                                            ind_cdiac) - 1 + 1,
                                                          i - 1]

# SHIFT
if (data_LULCC[:3] == 'LUH'):
    for b1 in range(len(bio)):
        for b2 in range(b1, len(bio)):
            if os.path.isfile(
                    'data/LandUse_' + data_LULCC + '/#DATA.LandUse_' + data_LULCC + '_' + mod_LSNKcover + '.1501-2015_114reg1.SHIFT_' +
                    bio[b1] + '2' + bio[b2] + '.csv'):
                TMP = np.array([line for line in csv.reader(open(
                    'data/LandUse_LUH1/#DATA.LandUse_' + data_LULCC + '_' + mod_LSNKcover + '.1501-2015_114reg1.SHIFT_' +
                    bio[b1] + '2' + bio[b2] + '.csv', 'r'))], dtype=dty)
                for i in range(1, 114 + 1):
                    SHIFT[1:min(ind_final, ind_cdiac) + 1, regionJ_index[i], 0, kLUC,
                    regionI_index[i], biome_index[bio[b1]], biome_index[bio[b2]]] += TMP[
                                                                                     200:200 + min(
                                                                                         ind_final,
                                                                                         ind_cdiac) - 1 + 1,
                                                                                     i - 1]

# ========
# 2.2. RCP
# ========

# initialization of projected drivers
LUCproj = np.zeros(
    [ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_biome, nb_biome],
    dtype=dty)
HARVproj = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_biome],
                    dtype=dty)
SHIFTproj = np.zeros(
    [ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_biome, nb_biome],
    dtype=dty)

# projection of land-use data under RCP scenarios
# from [Hurtt et al., 2011] and [Meinshausen et al., 2011]

# LUC
if (scen_LULCC[:3] == 'RCP') & (ind_final > ind_cdiac):
    for b1 in range(len(bio)):
        for b2 in range(len(bio)):
            if os.path.isfile(
                    'data/LandUse_RCP/#DATA.LandUse_RCP_' + mod_LSNKcover + '.2006-2100_114reg1.rcp' +
                    scen_LULCC[3] + scen_LULCC[5] + '_LUC_' + bio[b1] + '2' + bio[
                        b2] + '.csv'):
                TMP = np.array([line for line in csv.reader(open(
                    'data/LandUse_RCP/#DATA.LandUse_RCP_' + mod_LSNKcover + '.2006-2100_114reg1.rcp' +
                    scen_LULCC[3] + scen_LULCC[5] + '_LUC_' + bio[b1] + '2' + bio[
                        b2] + '.csv', 'r'))], dtype=dty)
                for i in range(1, 114 + 1):
                    LUCproj[306:min(ind_final, 400) + 1, regionJ_index[i], 0, kLUC,
                    regionI_index[i], biome_index[bio[b1]], biome_index[bio[b2]]] += TMP[
                                                                                     :min(
                                                                                         ind_final,
                                                                                         400) - 306 + 1,
                                                                                     i - 1]

# HARV
if (scen_LULCC[:3] == 'RCP') & (ind_final > ind_cdiac):
    for b in range(len(bio)):
        if os.path.isfile(
                'data/LandUse_RCP/#DATA.LandUse_RCP_' + mod_LSNKcover + '.2006-2100_114reg1.rcp' +
                scen_LULCC[3] + scen_LULCC[5] + '_HARV_' + bio[b] + '.csv'):
            TMP = np.array([line for line in csv.reader(open(
                'data/LandUse_RCP/#DATA.LandUse_RCP_' + mod_LSNKcover + '.2006-2100_114reg1.rcp' +
                scen_LULCC[3] + scen_LULCC[5] + '_HARV_' + bio[b] + '.csv', 'r'))],
                           dtype=dty)
            for i in range(1, 114 + 1):
                HARVproj[306:min(ind_final, 400) + 1, regionJ_index[i], 0, kLUC,
                regionI_index[i], biome_index[bio[b]]] += TMP[
                                                          :min(ind_final, 400) - 306 + 1,
                                                          i - 1]

# SHIFT
if (scen_LULCC[:3] == 'RCP') & (ind_final > ind_cdiac):
    for b1 in range(len(bio)):
        for b2 in range(b1, len(bio)):
            if os.path.isfile(
                    'data/LandUse_RCP/#DATA.LandUse_RCP_' + mod_LSNKcover + '.2006-2100_114reg1.rcp' +
                    scen_LULCC[3] + scen_LULCC[5] + '_SHIFT_' + bio[b1] + '2' + bio[
                        b2] + '.csv'):
                TMP = np.array([line for line in csv.reader(open(
                    'data/LandUse_RCP/#DATA.LandUse_RCP_' + mod_LSNKcover + '.2006-2100_114reg1.rcp' +
                    scen_LULCC[3] + scen_LULCC[5] + '_SHIFT_' + bio[b1] + '2' + bio[
                        b2] + '.csv', 'r'))], dtype=dty)
                for i in range(1, 114 + 1):
                    SHIFTproj[306:min(ind_final, 400) + 1, regionJ_index[i], 0, kLUC,
                    regionI_index[i], biome_index[bio[b1]], biome_index[bio[b2]]] += TMP[
                                                                                     :min(
                                                                                         ind_final,
                                                                                         400) - 306 + 1,
                                                                                     i - 1]

                # ==================
# 2.A. FINAL DATASET
# ==================

# datasets mixed following various criteria
for VAR in ['LUC', 'HARV', 'SHIFT']:

    # stop emissions
    if (scen_LULCC == 'stop') & (ind_final > ind_cdiac):
        exec(VAR + '[ind_cdiac+1:,...] = 0')

    # constant emissions
    elif (scen_LULCC == 'cst') & (ind_final > ind_cdiac):
        exec(VAR + '[ind_cdiac+1:,...] = ' + VAR + '[ind_cdiac,...][np.newaxis,...]')

        # RCP scenarios
    # always raw discontinuity
    elif (scen_LULCC[:3] == 'RCP') & (ind_final > ind_cdiac):
        exec(VAR + '[ind_cdiac+1:,...] = ' + VAR + 'proj[ind_cdiac+1:,...]')

# delete individual datasets
for VAR in ['LUC', 'HARV', 'SHIFT']:
    exec('del ' + VAR + 'proj')

##################################################
#   3. HALOGENATED COMPOUNDS
##################################################

EHFC = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_HFC],
                dtype=dty)  # {kt/yr}
EPFC = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_PFC],
                dtype=dty)  # {kt/yr}
EODS = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_ODS],
                dtype=dty)  # {kt/yr}

# ==========
# 3.1. EDGAR
# ==========

# load emissions from EDGAR v4.2
# see [JRC, 2011]

# HFCs
EHFCedgar = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_HFC],
                     dtype=dty)
for VAR in HFC:
    TMP = np.array([line for line in csv.reader(open(
        'data/EHaloComp_EDGAR/#DATA.EHaloComp_EDGAR.1970-' + str(
            1700 + ind_edgar) + '_114reg0.E' + VAR + '.csv', 'r'))], dtype=dty)
    for i in range(114 + 1):
        EHFCedgar[270:ind_edgar + 1, regionJ_index[i], 0, kindGHG_index['HFC'],
        regionI_index[i], HFC.index(VAR)] += TMP[:ind_edgar - 270 + 1, i]

# PFCs
EPFCedgar = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_PFC],
                     dtype=dty)
for VAR in PFC:
    TMP = np.array([line for line in csv.reader(open(
        'data/EHaloComp_EDGAR/#DATA.EHaloComp_EDGAR.1970-' + str(
            1700 + ind_edgar) + '_114reg0.E' + VAR + '.csv', 'r'))], dtype=dty)
    for i in range(114 + 1):
        EPFCedgar[270:ind_edgar + 1, regionJ_index[i], 0, kindGHG_index['PFC'],
        regionI_index[i], PFC.index(VAR)] += TMP[:ind_edgar - 270 + 1, i]

# =============
# 3.2. EDGAR-FT
# =============

# load emissions from EDGAR-FT v4.2-FT2010
# see [JRC, 2013]

# HFCs
EHFCeft = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_HFC],
                   dtype=dty)
for VAR in HFC:
    TMP = np.array([line for line in csv.reader(open(
        'data/EHaloComp_EDGAR-FT/#DATA.EHaloComp_EDGAR-FT.2008-2010_114reg0.E' + VAR + '.csv',
        'r'))], dtype=dty)
    for i in range(114 + 1):
        EHFCeft[308:310 + 1, regionJ_index[i], 0, kindGHG_index['HFC'], regionI_index[i],
        HFC.index(VAR)] += TMP[:310 - 308 + 1, i]

# PFCs
EPFCeft = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_PFC],
                   dtype=dty)
for VAR in PFC:
    TMP = np.array([line for line in csv.reader(open(
        'data/EHaloComp_EDGAR-FT/#DATA.EHaloComp_EDGAR-FT.2008-2010_114reg0.E' + VAR + '.csv',
        'r'))], dtype=dty)
    for i in range(114 + 1):
        EPFCeft[308:310 + 1, regionJ_index[i], 0, kindGHG_index['PFC'], regionI_index[i],
        PFC.index(VAR)] += TMP[:310 - 308 + 1, i]

# ==========
# 3.3. CMIP5
# ==========

# ODSs
EODScmip5 = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI, nb_ODS],
                     dtype=dty)
TMP = np.array([line for line in csv.reader(
    open('data/EHaloComp_CMIP5/#DATA.EHaloComp_CMIP5.1765-2005_(16ghg).EODS.csv', 'r'))][
               1:], dtype=dty)
lgd = [line for line in csv.reader(
    open('data/EHaloComp_CMIP5/#DATA.EHaloComp_CMIP5.1765-2005_(16ghg).EODS.csv', 'r'))][
    0]
for x in range(len(lgd)):
    EODScmip5[65:305 + 1, 0, 0, kindGHG_index['ODS'], 0, ODS.index(lgd[x])] = TMP[:min(
        ind_final, 305) - 65 + 1, x]

# extend dataset following RCP unique projection
TMP = np.array([line for line in csv.reader(
    open('data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_(16ghg).rcp85_EODS.csv',
         'r'))][1:], dtype=dty)
lgd = [line for line in csv.reader(
    open('data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_(16ghg).rcp85_EODS.csv',
         'r'))][0]
for x in range(len(lgd)):
    EODScmip5[306:ind_cdiac + 1, 0, 0, kindGHG_index['ODS'], 0, ODS.index(lgd[x])] = TMP[
                                                                                     6:min(
                                                                                         ind_final,
                                                                                         ind_cdiac) - 300 + 1,
                                                                                     x]

# ========
# 3.4. RCP
# ========

# initialization of projected drivers
for VAR in ['HFC', 'PFC', 'ODS']:
    exec(
        'E' + VAR + 'proj = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI,nb_' + VAR + '], dtype=dty)')

# projection of emissions under RCP scenarios
# from [Meinshausen et al., 2011]

# HFCs
if (scen_Ehalo[:3] == 'RCP') & (ind_final > ind_cdiac):
    for VAR in HFC:
        if (os.path.isfile(
                'data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_5reg0.rcp' + scen_Ehalo[
                    3] + scen_Ehalo[5] + '_E' + VAR + '.csv')):
            TMP = np.array([line for line in csv.reader(open(
                'data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_5reg0.rcp' + scen_Ehalo[
                    3] + scen_Ehalo[5] + '_E' + VAR + '.csv', 'r'))], dtype=dty)
            for i in range(4 + 1):
                if (mod_regionI == 'RCP5') & (mod_regionJ == 'RCP5'):
                    EHFCproj[300:min(ind_final, 400) + 1, i, 0, kindGHG_index['HFC'], i,
                    HFC.index(VAR)] += TMP[:min(ind_final, 400) - 300 + 1, i]
                elif (mod_regionI == 'RCP5') & (mod_regionJ != 'RCP5'):
                    EHFCproj[300:min(ind_final, 400) + 1, 0, 0, kindGHG_index['HFC'], i,
                    HFC.index(VAR)] += TMP[:min(ind_final, 400) - 300 + 1, i]
                elif (mod_regionI != 'RCP5') & (mod_regionJ == 'RCP5'):
                    EHFCproj[300:min(ind_final, 400) + 1, i, 0, kindGHG_index['HFC'], 0,
                    HFC.index(VAR)] += TMP[:min(ind_final, 400) - 300 + 1, i]
                else:
                    EHFCproj[300:min(ind_final, 400) + 1, 0, 0, kindGHG_index['HFC'], 0,
                    HFC.index(VAR)] += TMP[:min(ind_final, 400) - 300 + 1, i]

# PFCs
if (scen_Ehalo[:3] == 'RCP') & (ind_final > ind_cdiac):
    for VAR in PFC:
        if (os.path.isfile(
                'data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_5reg0.rcp' + scen_Ehalo[
                    3] + scen_Ehalo[5] + '_E' + VAR + '.csv')):
            TMP = np.array([line for line in csv.reader(open(
                'data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_5reg0.rcp' + scen_Ehalo[
                    3] + scen_Ehalo[5] + '_E' + VAR + '.csv', 'r'))], dtype=dty)
            for i in range(4 + 1):
                if (mod_regionI == 'RCP5') & (mod_regionJ == 'RCP5'):
                    EPFCproj[300:min(ind_final, 400) + 1, i, 0, kindGHG_index['PFC'], i,
                    PFC.index(VAR)] += TMP[:min(ind_final, 400) - 300 + 1, i]
                elif (mod_regionI == 'RCP5') & (mod_regionJ != 'RCP5'):
                    EPFCproj[300:min(ind_final, 400) + 1, 0, 0, kindGHG_index['PFC'], i,
                    PFC.index(VAR)] += TMP[:min(ind_final, 400) - 300 + 1, i]
                elif (mod_regionI != 'RCP5') & (mod_regionJ == 'RCP5'):
                    EPFCproj[300:min(ind_final, 400) + 1, i, 0, kindGHG_index['PFC'], 0,
                    PFC.index(VAR)] += TMP[:min(ind_final, 400) - 300 + 1, i]
                else:
                    EPFCproj[300:min(ind_final, 400) + 1, 0, 0, kindGHG_index['PFC'], 0,
                    PFC.index(VAR)] += TMP[:min(ind_final, 400) - 300 + 1, i]

# ODSs
if (scen_Ehalo[:3] == 'RCP') & (ind_final > ind_cdiac):
    TMP = np.array([line for line in csv.reader(open(
        'data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_(16ghg).rcp' + scen_Ehalo[3] +
        scen_Ehalo[5] + '_EODS.csv', 'r'))][1:], dtype=dty)
    lgd = [line for line in csv.reader(open(
        'data/EHaloComp_RCP/#DATA.EHaloComp_RCP.2000-2100_(16ghg).rcp' + scen_Ehalo[3] +
        scen_Ehalo[5] + '_EODS.csv', 'r'))][0]
    for x in range(len(lgd)):
        EODSproj[300:min(ind_final, 400) + 1, 0, 0, kindGHG_index['ODS'], 0,
        ODS.index(lgd[x])] = TMP[:min(ind_final, 400) - 300 + 1, x]

# =================
# 3.A. PAST DATASET
# =================

# datasets mixed following trends
for VAR in ['HFC', 'PFC', 'ODS']:

    if VAR in ['HFC']:

        # with EDGAR as reference
        if (data_Ehalo == 'EDGAR'):
            exec('E' + VAR + 'past = E' + VAR + 'edgar.copy()')
            # follow EDGAR-FT variations after 2008
            for t in range(ind_edgar + 1, ind_cdiac + 1):
                exec(
                    'E' + VAR + 'past[t,...] = E' + VAR + 'past[t-1,...] * E' + VAR + 'eft[t,...]/E' + VAR + 'eft[t-1,...]')
                exec(
                    'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
                exec(
                    'E' + VAR + 'past[t,...] *= np.sum(E' + VAR + 'past[t-1,...])/np.sum(E' + VAR + 'past[t,...]) * np.sum(E' + VAR + 'eft[t,...])/np.sum(E' + VAR + 'eft[t-1,...])')
            exec(
                'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
            # quadratic extrapolation before 1970
            # starting year of emission based on [Meinshausen et al., 2011]
            for x in range(nb_HFC):
                if (HFC[x] == 'HFC23'):
                    for t in range(230, 270):
                        exec(
                            'E' + VAR + 'past[t,...,x] = E' + VAR + 'past[270,...,x] * ((t-230)/(270-230.))**2')

    elif VAR in ['PFC']:

        # with EDGAR as reference
        if (data_Ehalo == 'EDGAR'):
            exec('E' + VAR + 'past = E' + VAR + 'edgar.copy()')
            # follow EDGAR-FT variations after 2008
            for t in range(ind_edgar + 1, ind_cdiac + 1):
                exec(
                    'E' + VAR + 'past[t,...] = E' + VAR + 'past[t-1,...] * E' + VAR + 'eft[t,...]/E' + VAR + 'eft[t-1,...]')
                exec(
                    'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
                exec(
                    'E' + VAR + 'past[t,...] *= np.sum(E' + VAR + 'past[t-1,...])/np.sum(E' + VAR + 'past[t,...]) * np.sum(E' + VAR + 'eft[t,...])/np.sum(E' + VAR + 'eft[t-1,...])')
            exec(
                'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
            # quadratic extrapolation before 1970
            # starting year of emission based on [Meinshausen et al., 2011]
            for x in range(nb_PFC):
                if (PFC[x] == 'SF6'):
                    for t in range(250, 270):
                        exec(
                            'E' + VAR + 'past[t,...,x] = E' + VAR + 'past[270,...,x] * ((t-250.)/(270-250.))**2')
                elif (PFC[x] == 'CF4'):
                    for t in range(222, 270):
                        exec(
                            'E' + VAR + 'past[t,...,x] = E' + VAR + 'past[270,...,x] * ((t-222.)/(270-222.))**2')
                elif (PFC[x] == 'C2F6'):
                    for t in range(189, 270):
                        exec(
                            'E' + VAR + 'past[t,...,x] = E' + VAR + 'past[270,...,x] * ((t-189)/(270-189.))**2')

    elif VAR in ['ODS']:

        # with CMIP5 as reference
        if (data_Ehalo == data_Ehalo):
            exec('E' + VAR + 'past = E' + VAR + 'cmip5.copy()')
            # linear extrapolation before 1765
            for t in range(50, 65):
                exec(
                    'E' + VAR + 'past[t,...] = E' + VAR + 'past[65,...] * (t-50)/float(65-50)')

    # cut past dataset to right length
    exec(
        'E' + VAR + '[:min(ind_cdiac,ind_final)+1,...] = E' + VAR + 'past[:min(ind_cdiac,ind_final)+1,...]')

# ==================
# 3.B. FINAL DATASET
# ==================

# datasets mixed following various criteria
for VAR in ['HFC', 'PFC', 'ODS']:

    # stop emissions
    if (scen_Ehalo == 'stop') & (ind_final > ind_cdiac):
        exec('E' + VAR + '[ind_cdiac+1:,...] = 0')

    # constant emissions
    elif (scen_Ehalo == 'cst') & (ind_final > ind_cdiac):
        exec(
            'E' + VAR + '[ind_cdiac+1:,...] = E' + VAR + '[ind_cdiac,...][np.newaxis,...]')

        # RCP scenarios
    elif (scen_Ehalo[:3] == 'RCP') & (ind_final > ind_cdiac):

        # raw discontinuity
        if (mod_DATAscen == 'raw'):
            exec('E' + VAR + '[ind_cdiac+1:,...] = E' + VAR + 'proj[ind_cdiac+1:,...]')

        # offset at transition point
        elif (mod_DATAscen == 'offset'):
            exec(
                'E' + VAR + '[ind_cdiac+1:,...] = E' + VAR + 'proj[ind_cdiac+1:,...] - E' + VAR + 'proj[ind_cdiac,...] + E' + VAR + '[ind_cdiac,...]')
            for t in range(ind_cdiac + 1, ind_final + 1):
                exec('def_regI = bool(np.sum(E' + VAR + 'proj[t,:,...,1:]))')
                exec('def_regJ = bool(np.sum(E' + VAR + 'proj[t,1:,...,:]))')
                if not def_regI:
                    exec('E' + VAR + '[t,:,...,0] += np.sum(E' + VAR + '[t,:,...,1:],-1)')
                    exec('E' + VAR + '[t,:,...,1:] = 0')
                if not def_regJ:
                    exec('E' + VAR + '[t,0,...,:] += np.sum(E' + VAR + '[t,1:,...,:],0)')
                    exec('E' + VAR + '[t,1:,...,:] = 0')

                    # linear transition over N years
        elif (mod_DATAscen[:6] == 'smooth'):
            N = int(mod_DATAscen[6:])
            if (ind_final >= ind_cdiac + N):
                for t in range(ind_cdiac + 1, ind_cdiac + N):
                    exec(
                        'E' + VAR + '[t,...] = (1-(t-ind_cdiac)/float(N)) * E' + VAR + '[ind_cdiac,...] + (t-ind_cdiac)/float(N) * E' + VAR + 'proj[ind_cdiac+N,...]')
                    exec('def_regI = bool(np.sum(E' + VAR + 'proj[t,:,...,1:]))')
                    exec('def_regJ = bool(np.sum(E' + VAR + 'proj[t,1:,...,:]))')
                    if not def_regI:
                        exec(
                            'E' + VAR + '[t,:,...,0] += np.sum(E' + VAR + '[t,:,...,1:],-1)')
                        exec('E' + VAR + '[t,:,...,1:] = 0')
                    if not def_regJ:
                        exec(
                            'E' + VAR + '[t,0,...,:] += np.sum(E' + VAR + '[t,1:,...,:],0)')
                        exec('E' + VAR + '[t,1:,...,:] = 0')
                exec(
                    'E' + VAR + '[ind_cdiac+N:,...] = E' + VAR + 'proj[ind_cdiac+N:,...]')

        # follow trends of projection
        elif (mod_DATAscen == 'trends'):
            for t in range(ind_cdiac + 1, ind_final + 1):
                exec('def_regI = bool(np.sum(E' + VAR + 'proj[t,:,...,1:]))')
                exec('def_regJ = bool(np.sum(E' + VAR + 'proj[t,1:,...,:]))')
                if def_regI and def_regJ:
                    exec(
                        'E' + VAR + '[t,...] = E' + VAR + '[t-1,...] * E' + VAR + 'proj[t,...]/E' + VAR + 'proj[t-1,...]')
                    exec(
                        'E' + VAR + '[np.isnan(E' + VAR + ')|np.isinf(E' + VAR + ')] = 0')
                    exec(
                        'E' + VAR + '[t,...] *= np.sum(E' + VAR + '[t-1,...])/np.sum(E' + VAR + '[t,...]) * np.sum(E' + VAR + 'proj[t,...])/np.sum(E' + VAR + 'proj[t-1,...])')
                elif not def_regI and def_regJ:
                    exec(
                        'E' + VAR + '[t,:,...,0] = np.sum(E' + VAR + '[t-1,:,...,:],0) * np.sum(E' + VAR + 'proj[t,:,...,:],-1)/np.sum(E' + VAR + 'proj[t-1,:,...,:],-1)')
                    exec(
                        'E' + VAR + '[np.isnan(E' + VAR + ')|np.isinf(E' + VAR + ')] = 0')
                    exec(
                        'E' + VAR + '[t,...] *= np.sum(E' + VAR + '[t-1,...])/np.sum(E' + VAR + '[t,...]) * np.sum(E' + VAR + 'proj[t,...])/np.sum(E' + VAR + 'proj[t-1,...])')
                elif def_regI and not def_regJ:
                    exec(
                        'E' + VAR + '[t,0,...,:] = np.sum(E' + VAR + '[t-1,:,...,:],0) * np.sum(E' + VAR + 'proj[t,:,...,:],0)/np.sum(E' + VAR + 'proj[t-1,:,...,:],0)')
                    exec(
                        'E' + VAR + '[np.isnan(E' + VAR + ')|np.isinf(E' + VAR + ')] = 0')
                    exec(
                        'E' + VAR + '[t,...] *= np.sum(E' + VAR + '[t-1,...])/np.sum(E' + VAR + '[t,...]) * np.sum(E' + VAR + 'proj[t,...])/np.sum(E' + VAR + 'proj[t-1,...])')
                elif not def_regI and not def_regJ:
                    exec(
                        'E' + VAR + '[t,0,...,0] = np.sum(np.sum(E' + VAR + '[t-1,:,...,:],-1),0) * np.sum(np.sum(E' + VAR + 'proj[t,:,...,:],-1),0)/np.sum(np.sum(E' + VAR + 'proj[t-1,:,...,:],-1),0)')
            exec('E' + VAR + '[np.isnan(E' + VAR + ')|np.isinf(E' + VAR + ')] = 0')

# delete individual datasets
for VAR in ['HFC', 'PFC']:
    exec('del E' + VAR + 'edgar,E' + VAR + 'eft,E' + VAR + 'past,E' + VAR + 'proj')
for VAR in ['ODS']:
    exec('del E' + VAR + 'cmip5,E' + VAR + 'past,E' + VAR + 'proj')

##################################################
#   4. SHORT-LIVED SPECIES
##################################################

ENOX = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                dtype=dty)  # {TgN/yr}
ECO = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
               dtype=dty)  # {TgC/yr}
EVOC = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                dtype=dty)  # {Tg/yr}
ESO2 = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                dtype=dty)  # {TgS/yr}
ENH3 = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                dtype=dty)  # {TgN/yr}
EOC = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
               dtype=dty)  # {Tg/yr}
EBC = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
               dtype=dty)  # {Tg/yr}

# ==========
# 4.1. EDGAR
# ==========

# load emissions from EDGAR v4.2
# see [JRC, 2011]

# OzoPrec
for VAR in ['NOX', 'CO', 'VOC']:
    exec(
        'E' + VAR + 'edgar = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)')
    for s in range(1, len(sec_accmip) - 2):
        if os.path.isfile('data/EOzoPrec_EDGAR/#DATA.EOzoPrec_EDGAR.1970-' + str(
                1700 + ind_edgar) + '_114reg0.E' + VAR + '_' + sec_accmip[s] + '.csv'):
            TMP = np.array([line for line in csv.reader(open(
                'data/EOzoPrec_EDGAR/#DATA.EOzoPrec_EDGAR.1970-' + str(
                    1700 + ind_edgar) + '_114reg0.E' + VAR + '_' + sec_accmip[s] + '.csv',
                'r'))], dtype=dty)
            for i in range(114 + 1):
                exec(
                    'E' + VAR + 'edgar[270:ind_edgar+1,regionJ_index[i],0,kindCHI_index[VAR],regionI_index[i]] += TMP[:ind_edgar-270+1,i]')

# AeroPrec
for VAR in ['SO2', 'NH3']:
    exec(
        'E' + VAR + 'edgar = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)')
    for s in range(1, len(sec_accmip) - 2):
        if os.path.isfile('data/EAeroPrec_EDGAR/#DATA.EAeroPrec_EDGAR.1970-' + str(
                1700 + ind_edgar) + '_114reg0.E' + VAR + '_' + sec_accmip[s] + '.csv'):
            TMP = np.array([line for line in csv.reader(open(
                'data/EAeroPrec_EDGAR/#DATA.EAeroPrec_EDGAR.1970-' + str(
                    1700 + ind_edgar) + '_114reg0.E' + VAR + '_' + sec_accmip[s] + '.csv',
                'r'))], dtype=dty)
            for i in range(114 + 1):
                exec(
                    'E' + VAR + 'edgar[270:ind_edgar+1,regionJ_index[i],0,kindAER_index[VAR],regionI_index[i]] += TMP[:ind_edgar-270+1,i]')

# PM10 (proxy for OC/BC)
EPM10edgar = np.zeros([ind_cdiac + 1, nb_regionJ, nb_sector, nb_kind, nb_regionI],
                      dtype=dty)
for s in range(1, len(sec_accmip) - 2):
    if os.path.isfile('data/EAeroPrec_EDGAR/#DATA.EAeroPrec_EDGAR.1970-' + str(
            1700 + ind_edgar) + '_114reg0.EPM10_' + sec_accmip[s] + '.csv'):
        TMP = np.array([line for line in csv.reader(open(
            'data/EAeroPrec_EDGAR/#DATA.EAeroPrec_EDGAR.1970-' + str(
                1700 + ind_edgar) + '_114reg0.EPM10_' + sec_accmip[s] + '.csv', 'r'))],
                       dtype=dty)
        for i in range(114 + 1):
            EPM10edgar[270:ind_edgar + 1, regionJ_index[i], 0, kindAER_index['OC'],
            regionI_index[i]] += TMP[:ind_edgar - 270 + 1, i] / 2
            EPM10edgar[270:ind_edgar + 1, regionJ_index[i], 0, kindAER_index['BC'],
            regionI_index[i]] += TMP[:ind_edgar - 270 + 1, i] / 2

# ===============
# 4.2. EDGAR-HTAP
# ===============

# load emissions from EDGAR-HTAP v2
# see [JRC, 2013]

# OzoPrec
for VAR in ['NOX', 'CO', 'VOC']:
    exec(
        'E' + VAR + 'ehtap = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)')
    for s in range(1, len(sec_accmip) - 2):
        if os.path.isfile(
                'data/EOzoPrec_EDGAR-HTAP/#DATA.EOzoPrec_EDGAR-HTAP.2008-2010_114reg0.E' + VAR + '_' +
                sec_accmip[s] + '.csv'):
            TMP = np.array([line for line in csv.reader(open(
                'data/EOzoPrec_EDGAR-HTAP/#DATA.EOzoPrec_EDGAR-HTAP.2008-2010_114reg0.E' + VAR + '_' +
                sec_accmip[s] + '.csv', 'r'))], dtype=dty)
            for i in range(114 + 1):
                exec(
                    'E' + VAR + 'ehtap[308:310+1,regionJ_index[i],0,kindCHI_index[VAR],regionI_index[i]] += TMP[:ind_edgar-270+1,i]')

# AeroPrec
for VAR in ['SO2', 'NH3', 'OC', 'BC']:
    exec(
        'E' + VAR + 'ehtap = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)')
    for s in range(1, len(sec_accmip) - 2):
        if os.path.isfile(
                'data/EAeroPrec_EDGAR-HTAP/#DATA.EAeroPrec_EDGAR-HTAP.2008-2010_114reg0.E' + VAR + '_' +
                sec_accmip[s] + '.csv'):
            TMP = np.array([line for line in csv.reader(open(
                'data/EAeroPrec_EDGAR-HTAP/#DATA.EAeroPrec_EDGAR-HTAP.2008-2010_114reg0.E' + VAR + '_' +
                sec_accmip[s] + '.csv', 'r'))], dtype=dty)
            for i in range(114 + 1):
                exec(
                    'E' + VAR + 'ehtap[308:310+1,regionJ_index[i],0,kindAER_index[VAR],regionI_index[i]] += TMP[:ind_edgar-270+1,i]')

# ===========
# 4.3. ACCMIP
# ===========

# load emissions from ACCMIP
# from [Lamarque et al., 2010]

# OzoPrec
for VAR in ['NOX', 'CO', 'VOC']:
    exec(
        'E' + VAR + 'accmip = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)')
    exec(
        'p_E' + VAR + '_bio = np.zeros([nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)')
    for s in range(1, len(sec_accmip) - 2):
        if os.path.isfile(
                'data/EOzoPrec_ACCMIP/#DATA.EOzoPrec_ACCMIP.1850-2000_114reg0.E' + VAR + '_' +
                sec_accmip[s] + '.csv'):
            TMP = np.array([line for line in csv.reader(open(
                'data/EOzoPrec_ACCMIP/#DATA.EOzoPrec_ACCMIP.1850-2000_114reg0.E' + VAR + '_' +
                sec_accmip[s] + '.csv', 'r'))], dtype=dty)
            for i in range(114 + 1):
                exec(
                    'E' + VAR + 'accmip[150:300+1,regionJ_index[i],0,kindCHI_index[VAR],regionI_index[i]] += TMP[:300-150+1,i]')
                if sec_accmip[s] in ['agr', 'awb', 'wst']:
                    exec(
                        'p_E' + VAR + '_bio[regionJ_index[i],0,kindCHI_index[VAR],regionI_index[i]] += TMP[0,i]')
    exec('p_E' + VAR + '_bio /= E' + VAR + 'accmip[150]')
    exec(
        'p_E' + VAR + '_bio[np.isnan(p_E' + VAR + '_bio)|np.isinf(p_E' + VAR + '_bio)] = 0')

# AeroPrec
for VAR in ['SO2', 'NH3', 'OC', 'BC']:
    exec(
        'E' + VAR + 'accmip = np.zeros([ind_cdiac+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)')
    exec(
        'p_E' + VAR + '_bio = np.zeros([nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)')
    for s in range(1, len(sec_accmip) - 2):
        if os.path.isfile(
                'data/EAeroPrec_ACCMIP/#DATA.EAeroPrec_ACCMIP.1850-2000_114reg0.E' + VAR + '_' +
                sec_accmip[s] + '.csv'):
            TMP = np.array([line for line in csv.reader(open(
                'data/EAeroPrec_ACCMIP/#DATA.EAeroPrec_ACCMIP.1850-2000_114reg0.E' + VAR + '_' +
                sec_accmip[s] + '.csv', 'r'))], dtype=dty)
            for i in range(114 + 1):
                exec(
                    'E' + VAR + 'accmip[150:300+1,regionJ_index[i],0,kindAER_index[VAR],regionI_index[i]] += TMP[:300-150+1,i]')
                if sec_accmip[s] in ['agr', 'awb', 'wst']:
                    exec(
                        'p_E' + VAR + '_bio[regionJ_index[i],0,kindAER_index[VAR],regionI_index[i]] += TMP[0,i]')
    exec('p_E' + VAR + '_bio /= E' + VAR + 'accmip[150]')
    exec(
        'p_E' + VAR + '_bio[np.isnan(p_E' + VAR + '_bio)|np.isinf(p_E' + VAR + '_bio)] = 0')

# =========
# 4.4. SRES
# =========

# initialization of projected drivers
for VAR in ['NOX', 'CO', 'VOC'] + ['SO2', 'NH3', 'OC', 'BC']:
    exec(
        'E' + VAR + 'proj = np.zeros([ind_final+1,nb_regionJ,nb_sector,nb_kind,nb_regionI], dtype=dty)')

# projection of emissions under SRES scenarios
# from [IPCC, 2000]

# OzoPrec
for VAR in ['NOX', 'CO', 'VOC']:
    exec('scen = scen_E' + VAR)
    if (scen[:4] == 'SRES') & (ind_final > ind_cdiac):
        TMP = np.array([line for line in csv.reader(open(
            'data/EOzoPrec_SRES/#DATA.EOzoPrec_SRES.2000-2100_4reg0.' + scen[
                                                                        5:] + '_E' + VAR + '.csv',
            'r'))], dtype=dty)
        for i in range(4 + 1):
            if (mod_regionI == 'SRES4') & (mod_regionJ == 'SRES4'):
                exec(
                    'E' + VAR + 'proj[300:min(ind_final,400)+1,i,0,kindCHI_index[VAR],i] += TMP[:min(ind_final,400)-300+1,i]')
            elif (mod_regionI == 'SRES4') & (mod_regionJ != 'SRES4'):
                exec(
                    'E' + VAR + 'proj[300:min(ind_final,400)+1,0,0,kindCHI_index[VAR],i] += TMP[:min(ind_final,400)-300+1,i]')
            elif (mod_regionI != 'SRES4') & (mod_regionJ == 'SRES4'):
                exec(
                    'E' + VAR + 'proj[300:min(ind_final,400)+1,i,0,kindCHI_index[VAR],0] += TMP[:min(ind_final,400)-300+1,i]')
            elif (mod_regionI != 'SRES4') & (mod_regionJ != 'SRES4'):
                exec(
                    'E' + VAR + 'proj[300:min(ind_final,400)+1,0,0,kindCHI_index[VAR],0] += TMP[:min(ind_final,400)-300+1,i]')

# AeroPrec (SO2 only)
if (scen_ESO2[:4] == 'SRES') & (ind_final > ind_cdiac):
    TMP = np.array([line for line in csv.reader(open(
        'data/EAeroPrec_SRES/#DATA.EAeroPrec_SRES.2000-2100_4reg0.' + scen_ESO2[
                                                                      5:] + '_ESO2.csv',
        'r'))], dtype=dty)
    for i in range(4 + 1):
        if (mod_regionI == 'SRES4') & (mod_regionJ == 'SRES4'):
            ESO2proj[300:min(ind_final, 400) + 1, i, 0, kindAER_index['SO2'], i] += TMP[
                                                                                    :min(
                                                                                        ind_final,
                                                                                        400) - 300 + 1,
                                                                                    i]
        elif (mod_regionI == 'SRES4') & (mod_regionJ != 'SRES4'):
            ESO2proj[300:min(ind_final, 400) + 1, 0, 0, kindAER_index['SO2'], i] += TMP[
                                                                                    :min(
                                                                                        ind_final,
                                                                                        400) - 300 + 1,
                                                                                    i]
        elif (mod_regionI != 'SRES4') & (mod_regionJ == 'SRES4'):
            ESO2proj[300:min(ind_final, 400) + 1, i, 0, kindAER_index['SO2'], 0] += TMP[
                                                                                    :min(
                                                                                        ind_final,
                                                                                        400) - 300 + 1,
                                                                                    i]
        elif (mod_regionI != 'SRES4') & (mod_regionJ != 'SRES4'):
            ESO2proj[300:min(ind_final, 400) + 1, 0, 0, kindAER_index['SO2'], 0] += TMP[
                                                                                    :min(
                                                                                        ind_final,
                                                                                        400) - 300 + 1,
                                                                                    i]

# ========
# 4.5. RCP
# ========

# projection of emissions under RCP scenarios
# from [Meinshausen et al., 2011]

# OzoPrec
for VAR in ['NOX', 'CO', 'VOC']:
    exec('scen = scen_E' + VAR)
    if (scen[:3] == 'RCP') & (ind_final > ind_cdiac):
        for s in range(1, len(sec_accmip) - 2):
            if os.path.isfile(
                    'data/EOzoPrec_RCP/#DATA.EOzoPrec_RCP.2000-2100_5reg0.rcp' + scen[3] +
                    scen[5] + '_E' + VAR + '_' + sec_accmip[s] + '.csv'):
                TMP = np.array([line for line in csv.reader(open(
                    'data/EOzoPrec_RCP/#DATA.EOzoPrec_RCP.2000-2100_5reg0.rcp' + scen[3] +
                    scen[5] + '_E' + VAR + '_' + sec_accmip[s] + '.csv', 'r'))],
                               dtype=dty)
                for i in range(5 + 1):
                    if (mod_regionI == 'RCP5') & (mod_regionJ == 'RCP5'):
                        exec(
                            'E' + VAR + 'proj[300:min(ind_final,400)+1,i,0,kindCHI_index[VAR],i] += TMP[:min(ind_final,400)-300+1,i]')
                    elif (mod_regionI == 'RCP5') & (mod_regionJ != 'RCP5'):
                        exec(
                            'E' + VAR + 'proj[300:min(ind_final,400)+1,0,0,kindCHI_index[VAR],i] += TMP[:min(ind_final,400)-300+1,i]')
                    elif (mod_regionI != 'RCP5') & (mod_regionJ == 'RCP5'):
                        exec(
                            'E' + VAR + 'proj[300:min(ind_final,400)+1,i,0,kindCHI_index[VAR],0] += TMP[:min(ind_final,400)-300+1,i]')
                    elif (mod_regionI != 'RCP5') & (mod_regionJ != 'RCP5'):
                        exec(
                            'E' + VAR + 'proj[300:min(ind_final,400)+1,0,0,kindCHI_index[VAR],0] += TMP[:min(ind_final,400)-300+1,i]')

# AeroPrec
for VAR in ['SO2', 'NH3', 'OC', 'BC']:
    exec('scen = scen_E' + VAR)
    if (scen[:3] == 'RCP') & (ind_final > ind_cdiac):
        for s in range(1, len(sec_accmip) - 2):
            if os.path.isfile(
                    'data/EAeroPrec_RCP/#DATA.EAeroPrec_RCP.2000-2100_5reg0.rcp' + scen[
                        3] + scen[5] + '_E' + VAR + '_' + sec_accmip[s] + '.csv'):
                TMP = np.array([line for line in csv.reader(open(
                    'data/EAeroPrec_RCP/#DATA.EAeroPrec_RCP.2000-2100_5reg0.rcp' + scen[
                        3] + scen[5] + '_E' + VAR + '_' + sec_accmip[s] + '.csv', 'r'))],
                               dtype=dty)
                for i in range(5 + 1):
                    if (mod_regionI == 'RCP5') & (mod_regionJ == 'RCP5'):
                        exec(
                            'E' + VAR + 'proj[300:min(ind_final,400)+1,i,0,kindAER_index[VAR],i] += TMP[:min(ind_final,400)-300+1,i]')
                    elif (mod_regionI == 'RCP5') & (mod_regionJ != 'RCP5'):
                        exec(
                            'E' + VAR + 'proj[300:min(ind_final,400)+1,0,0,kindAER_index[VAR],i] += TMP[:min(ind_final,400)-300+1,i]')
                    elif (mod_regionI != 'RCP5') & (mod_regionJ == 'RCP5'):
                        exec(
                            'E' + VAR + 'proj[300:min(ind_final,400)+1,i,0,kindAER_index[VAR],0] += TMP[:min(ind_final,400)-300+1,i]')
                    elif (mod_regionI != 'RCP5') & (mod_regionJ != 'RCP5'):
                        exec(
                            'E' + VAR + 'proj[300:min(ind_final,400)+1,0,0,kindAER_index[VAR],0] += TMP[:min(ind_final,400)-300+1,i]')

# =================
# 4.A. PAST DATASET
# =================

# datasets mixed following trends
for VAR in ['NOX', 'CO', 'VOC'] + ['SO2', 'NH3', 'OC', 'BC']:

    if VAR in ['NOX', 'CO', 'VOC'] + ['SO2', 'NH3']:
        exec('data = data_E' + VAR)

        # with EDGAR as reference
        if (data == 'EDGAR'):
            exec('E' + VAR + 'past = E' + VAR + 'edgar.copy()')
            # follow EDGAR-HTAP variations after 2008
            for t in range(ind_edgar + 1, ind_cdiac + 1):
                exec(
                    'E' + VAR + 'past[t,...] = E' + VAR + 'past[t-1,...] * E' + VAR + 'ehtap[t,...]/E' + VAR + 'ehtap[t-1,...]')
                exec(
                    'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
                exec(
                    'E' + VAR + 'past[t,...] *= np.sum(E' + VAR + 'past[t-1,...])/np.sum(E' + VAR + 'past[t,...]) * np.sum(E' + VAR + 'ehtap[t,...])/np.sum(E' + VAR + 'ehtap[t-1,...])')
            exec(
                'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
            # follow ACCMIP variations before 1970
            for t in range(150, 270)[::-1]:
                exec(
                    'E' + VAR + 'past[t,...] = E' + VAR + 'past[t+1,...] * E' + VAR + 'accmip[t,...]/E' + VAR + 'accmip[t+1,...]')
                exec(
                    'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
                exec(
                    'E' + VAR + 'past[t,...] *= np.sum(E' + VAR + 'past[t+1,...])/np.sum(E' + VAR + 'past[t,...]) * np.sum(E' + VAR + 'accmip[t,...])/np.sum(E' + VAR + 'accmip[t+1,...])')
            exec(
                'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
            # linear extrapolation before 1850
            for t in range(50, 150):
                exec(
                    'E' + VAR + 'past[t,...] = (1-p_E' + VAR + '_bio) * E' + VAR + 'past[150,...] * (t-50)/float(150-50)')
            for t in range(0, 150):
                exec(
                    'E' + VAR + 'past[t,...] += p_E' + VAR + '_bio * E' + VAR + 'past[150,...] * (t--200)/float(150--200)')
            exec('E' + VAR + '_0 = E' + VAR + 'past[50,...]')

        # with ACCMIP as reference
        elif (data == 'ACCMIP'):
            exec('E' + VAR + 'past = E' + VAR + 'accmip.copy()')
            # follow EDGAR variations after 2000
            for t in range(300 + 1, ind_edgar + 1):
                exec(
                    'E' + VAR + 'past[t,...] = E' + VAR + 'past[t-1,...] * E' + VAR + 'edgar[t,...]/E' + VAR + 'edgar[t-1,...]')
                exec(
                    'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
                exec(
                    'E' + VAR + 'past[t,...] *= np.sum(E' + VAR + 'past[t-1,...])/np.sum(E' + VAR + 'past[t,...]) * np.sum(E' + VAR + 'edgar[t,...])/np.sum(E' + VAR + 'edgar[t-1,...])')
            exec(
                'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
            # follow EDGAR-HTAP variations after 2008
            for t in range(ind_edgar + 1, ind_cdiac + 1):
                exec(
                    'E' + VAR + 'past[t,...] = E' + VAR + 'past[t-1,...] * E' + VAR + 'ehtap[t,...]/E' + VAR + 'ehtap[t-1,...]')
                exec(
                    'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
                exec(
                    'E' + VAR + 'past[t,...] *= np.sum(E' + VAR + 'past[t-1,...])/np.sum(E' + VAR + 'past[t,...]) * np.sum(E' + VAR + 'ehtap[t,...])/np.sum(E' + VAR + 'ehtap[t-1,...])')
            exec(
                'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
            # linear extrapolation before 1850
            for t in range(50, 150):
                exec(
                    'E' + VAR + 'past[t,...] = (1-p_E' + VAR + '_bio) * E' + VAR + 'past[150,...] * (t-50)/float(150-50)')
            for t in range(0, 150):
                exec(
                    'E' + VAR + 'past[t,...] += p_E' + VAR + '_bio * E' + VAR + 'past[150,...] * (t--200)/float(150--200)')
            exec('E' + VAR + '_0 = E' + VAR + 'past[50,...]')

    elif VAR in ['OC', 'BC']:
        exec('data = data_E' + VAR)

        # with ACCMIP as reference
        if (data == 'ACCMIP'):
            exec('E' + VAR + 'past = E' + VAR + 'accmip.copy()')
            # follow EDGAR (PM10) variations after 2000
            for t in range(300 + 1, ind_edgar + 1):
                exec(
                    'E' + VAR + 'past[t,...] = E' + VAR + 'past[t-1,...] * EPM10edgar[t,...]/EPM10edgar[t-1,...]')
                exec(
                    'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
                exec(
                    'E' + VAR + 'past[t,...] *= np.sum(E' + VAR + 'past[t-1,...])/np.sum(E' + VAR + 'past[t,...]) * np.sum(EPM10edgar[t,...])/np.sum(EPM10edgar[t-1,...])')
            exec(
                'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
            # follow EDGAR-HTAP variations after 2008
            for t in range(ind_edgar + 1, ind_cdiac + 1):
                exec(
                    'E' + VAR + 'past[t,...] = E' + VAR + 'past[t-1,...] * E' + VAR + 'ehtap[t,...]/E' + VAR + 'ehtap[t-1,...]')
                exec(
                    'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
                exec(
                    'E' + VAR + 'past[t,...] *= np.sum(E' + VAR + 'past[t-1,...])/np.sum(E' + VAR + 'past[t,...]) * np.sum(E' + VAR + 'ehtap[t,...])/np.sum(E' + VAR + 'ehtap[t-1,...])')
            exec(
                'E' + VAR + 'past[np.isnan(E' + VAR + 'past)|np.isinf(E' + VAR + 'past)] = 0')
            # linear extrapolation before 1850
            for t in range(50, 150):
                exec(
                    'E' + VAR + 'past[t,...] = (1-p_E' + VAR + '_bio) * E' + VAR + 'past[150,...] * (t-50)/float(150-50)')
            for t in range(0, 150):
                exec(
                    'E' + VAR + 'past[t,...] += p_E' + VAR + '_bio * E' + VAR + 'past[150,...] * (t--200)/float(150--200)')
            exec('E' + VAR + '_0 = E' + VAR + 'past[50,...]')

    # cut past dataset to right length
    exec(
        'E' + VAR + '[:min(ind_cdiac,ind_final)+1,...] = E' + VAR + 'past[:min(ind_cdiac,ind_final)+1,...]')

# ==================
# 4.B. FINAL DATASET
# ==================

# datasets mixed following various criteria
for VAR in ['NOX', 'CO', 'VOC'] + ['SO2', 'NH3', 'OC', 'BC']:
    exec('scen = scen_E' + VAR)

    # stop emissions
    if (scen == 'stop') & (ind_final > ind_cdiac):
        exec('E' + VAR + '[ind_cdiac+1:,...] = E' + VAR + '_0[np.newaxis,...]')

    # constant emissions
    elif (scen == 'cst') & (ind_final > ind_cdiac):
        exec(
            'E' + VAR + '[ind_cdiac+1:,...] = E' + VAR + '[ind_cdiac,...][np.newaxis,...]')

        # RCP or SRES scenarios
    elif ((scen[:4] == 'SRES') | (scen[:3] == 'RCP')) & (ind_final > ind_cdiac):

        # raw discontinuity
        if (mod_DATAscen == 'raw'):
            exec('E' + VAR + '[ind_cdiac+1:,...] = E' + VAR + 'proj[ind_cdiac+1:,...]')

        # offset at transition point
        elif (mod_DATAscen == 'offset'):
            exec(
                'E' + VAR + '[ind_cdiac+1:,...] = E' + VAR + 'proj[ind_cdiac+1:,...] - E' + VAR + 'proj[ind_cdiac,...] + E' + VAR + '[ind_cdiac,...]')
            for t in range(ind_cdiac + 1, ind_final + 1):
                exec('def_regI = bool(np.sum(E' + VAR + 'proj[t,:,...,1:]))')
                exec('def_regJ = bool(np.sum(E' + VAR + 'proj[t,1:,...,:]))')
                if not def_regI:
                    exec('E' + VAR + '[t,:,...,0] += np.sum(E' + VAR + '[t,:,...,1:],-1)')
                    exec('E' + VAR + '[t,:,...,1:] = 0')
                if not def_regJ:
                    exec('E' + VAR + '[t,0,...,:] += np.sum(E' + VAR + '[t,1:,...,:],0)')
                    exec('E' + VAR + '[t,1:,...,:] = 0')

                    # linear transition over N years
        elif (mod_DATAscen[:6] == 'smooth'):
            N = int(mod_DATAscen[6:])
            if (ind_final >= ind_cdiac + N):
                for t in range(ind_cdiac + 1, ind_cdiac + N):
                    exec(
                        'E' + VAR + '[t,...] = (1-(t-ind_cdiac)/float(N)) * E' + VAR + '[ind_cdiac,...] + (t-ind_cdiac)/float(N) * E' + VAR + 'proj[ind_cdiac+N,...]')
                    exec('def_regI = bool(np.sum(E' + VAR + 'proj[t,:,...,1:]))')
                    exec('def_regJ = bool(np.sum(E' + VAR + 'proj[t,1:,...,:]))')
                    if not def_regI:
                        exec(
                            'E' + VAR + '[t,:,...,0] += np.sum(E' + VAR + '[t,:,...,1:],-1)')
                        exec('E' + VAR + '[t,:,...,1:] = 0')
                    if not def_regJ:
                        exec(
                            'E' + VAR + '[t,0,...,:] += np.sum(E' + VAR + '[t,1:,...,:],0)')
                        exec('E' + VAR + '[t,1:,...,:] = 0')
                exec(
                    'E' + VAR + '[ind_cdiac+N:,...] = E' + VAR + 'proj[ind_cdiac+N:,...]')

        # follow trends of projection
        elif (mod_DATAscen == 'trends'):
            for t in range(ind_cdiac + 1, ind_final + 1):
                exec('def_regI = bool(np.sum(E' + VAR + 'proj[t,:,...,1:]))')
                exec('def_regJ = bool(np.sum(E' + VAR + 'proj[t,1:,...,:]))')
                if def_regI and def_regJ:
                    exec(
                        'E' + VAR + '[t,...] = E' + VAR + '[t-1,...] * E' + VAR + 'proj[t,...]/E' + VAR + 'proj[t-1,...]')
                    exec(
                        'E' + VAR + '[np.isnan(E' + VAR + ')|np.isinf(E' + VAR + ')] = 0')
                    exec(
                        'E' + VAR + '[t,...] *= np.sum(E' + VAR + '[t-1,...])/np.sum(E' + VAR + '[t,...]) * np.sum(E' + VAR + 'proj[t,...])/np.sum(E' + VAR + 'proj[t-1,...])')
                elif not def_regI and def_regJ:
                    exec(
                        'E' + VAR + '[t,:,...,0] = np.sum(E' + VAR + '[t-1,:,...,:],-1) * np.sum(E' + VAR + 'proj[t,:,...,:],-1)/np.sum(E' + VAR + 'proj[t-1,:,...,:],-1)')
                    exec(
                        'E' + VAR + '[np.isnan(E' + VAR + ')|np.isinf(E' + VAR + ')] = 0')
                    exec(
                        'E' + VAR + '[t,...] *= np.sum(E' + VAR + '[t-1,...])/np.sum(E' + VAR + '[t,...]) * np.sum(E' + VAR + 'proj[t,...])/np.sum(E' + VAR + 'proj[t-1,...])')
                elif def_regI and not def_regJ:
                    exec(
                        'E' + VAR + '[t,0,...,:] = np.sum(E' + VAR + '[t-1,:,...,:],0) * np.sum(E' + VAR + 'proj[t,:,...,:],0)/np.sum(E' + VAR + 'proj[t-1,:,...,:],0)')
                    exec(
                        'E' + VAR + '[np.isnan(E' + VAR + ')|np.isinf(E' + VAR + ')] = 0')
                    exec(
                        'E' + VAR + '[t,...] *= np.sum(E' + VAR + '[t-1,...])/np.sum(E' + VAR + '[t,...]) * np.sum(E' + VAR + 'proj[t,...])/np.sum(E' + VAR + 'proj[t-1,...])')
                elif not def_regI and not def_regJ:
                    exec(
                        'E' + VAR + '[t,0,...,0] = np.sum(np.sum(E' + VAR + '[t-1,:,...,:],-1),0) * np.sum(np.sum(E' + VAR + 'proj[t,:,...,:],-1),0)/np.sum(np.sum(E' + VAR + 'proj[t-1,:,...,:],-1),0)')
            exec('E' + VAR + '[np.isnan(E' + VAR + ')|np.isinf(E' + VAR + ')] = 0')

# delete individual datasets
for VAR in ['NOX', 'CO', 'VOC'] + ['SO2', 'NH3']:
    exec(
        'del E' + VAR + 'edgar,E' + VAR + 'ehtap,E' + VAR + 'accmip,E' + VAR + 'past,E' + VAR + 'proj,p_E' + VAR + '_bio')
for VAR in ['OC', 'BC']:
    exec(
        'del E' + VAR + 'ehtap,E' + VAR + 'accmip,E' + VAR + 'past,E' + VAR + 'proj,p_E' + VAR + '_bio')
for VAR in ['PM10']:
    exec('del E' + VAR + 'edgar')

##################################################
#   5. RADIATIVE FORCING
##################################################

RFcon = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind], dtype=dty)  # {W/m2}
RFvolc = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind], dtype=dty)  # {W/m2}
RFsolar = np.zeros([ind_final + 1, nb_regionJ, nb_sector, nb_kind], dtype=dty)  # {W/m2}

# ==================
# 5.1. Anthropogenic
# ==================

# load RF estimates from IPCC AR5
# from [IPCC, 2013] (annexe 2)
if (data_RFant == 'IPCC-AR5'):
    TMP = np.array([line for line in csv.reader(
        open('data/Historic_IPCC-AR5/#DATA.Historic_IPCC-AR5.1750-2011_(11for).RF.csv',
             'r'))][1:], dtype=dty)
    lgd = [line for line in csv.reader(
        open('data/Historic_IPCC-AR5/#DATA.Historic_IPCC-AR5.1750-2011_(11for).RF.csv',
             'r'))][0]
    # past estimates
    for x in range(len(lgd)):
        if (lgd[x] == 'Contrails'):
            RFcon[51:min(ind_final, ind_cdiac) + 1, 0, 0, kindRF_index['RFcon']] += TMP[
                                                                                    1:1 + min(
                                                                                        ind_final,
                                                                                        ind_cdiac) - 51 + 1,
                                                                                    x]
    # arbitrary extension
    if (scen_RFant == 'cst') & (ind_final > ind_cdiac):
        for x in range(len(lgd)):
            if (lgd[x] == 'Contrails'):
                RFcon[ind_cdiac + 1:, 0, 0, kindRF_index['RFcon']] = TMP[-1, x]

# ============
# 5.2. Natural
# ============

# load RF estimates from IPCC AR5
# from [IPCC, 2013] (annexe 2)
if (data_RFnat == 'IPCC-AR5'):
    TMP = np.array([line for line in csv.reader(
        open('data/Historic_IPCC-AR5/#DATA.Historic_IPCC-AR5.1750-2011_(11for).RF.csv',
             'r'))][1:], dtype=dty)
    lgd = [line for line in csv.reader(
        open('data/Historic_IPCC-AR5/#DATA.Historic_IPCC-AR5.1750-2011_(11for).RF.csv',
             'r'))][0]
    # past estimates
    for x in range(len(lgd)):
        if (lgd[x] == 'Volcano'):
            RFvolc[51:min(ind_final, ind_cdiac) + 1, 0, 0, kindRF_index['RFvol']] += TMP[
                                                                                     1:1 + min(
                                                                                         ind_final,
                                                                                         ind_cdiac) - 51 + 1,
                                                                                     x] - np.mean(
                TMP[1:, x])
        elif (lgd[x] == 'Solar'):
            RFsolar[51:min(ind_final, ind_cdiac) + 1, 0, 0, kindRF_index['RFsol']] += TMP[
                                                                                      1:1 + min(
                                                                                          ind_final,
                                                                                          ind_cdiac) - 51 + 1,
                                                                                      x]
    # arbitrary extension
    if (scen_RFnat == 'cst') & (ind_final > ind_cdiac):
        for x in range(len(lgd)):
            if (lgd[x] == 'Volcano'):
                RFvolc[ind_cdiac + 1:, 0, 0, kindRF_index['RFvol']] = 0.
            elif (lgd[x] == 'Solar'):
                RFsolar[ind_cdiac + 1:, 0, 0, kindRF_index['RFsol']] = np.mean(
                    TMP[-11:, x])

##################################################
#   B. FINAL DRIVERS
##################################################

# ==================
# B.1. PREINDUSTRIAL
# ==================

# reference emissions for preindustrial
for VAR in ['ECH4', 'EN2O'] + ['ENOX', 'ECO', 'EVOC', 'ESO2', 'ENH3', 'EOC', 'EBC']:
    exec(VAR + '[:,...] -= ' + VAR + '_0[np.newaxis,...]')

# force true 1750 preindustrial (drivers)
if PI_1750:
    for VAR in ['EFF', 'ECH4', 'EN2O'] + ['LUC', 'HARV', 'SHIFT'] + ['EHFC', 'EPFC',
                                                                     'EODS'] + ['ENOX',
                                                                                'ECO',
                                                                                'EVOC',
                                                                                'ESO2',
                                                                                'ENH3',
                                                                                'EOC',
                                                                                'EBC'] + [
                   'RFcon', 'RFvolc', 'RFsolar']:
        exec(VAR + '[:50+1] *= 0')

# ================
# B.2. ATTRIBUTION
# ================

# starting year of regional attribution
for VAR in ['EFF', 'ECH4', 'EN2O'] + ['LUC', 'HARV', 'SHIFT'] + ['EHFC', 'EPFC',
                                                                 'EODS'] + ['ENOX', 'ECO',
                                                                            'EVOC',
                                                                            'ESO2',
                                                                            'ENH3', 'EOC',
                                                                            'EBC'] + [
               'RFcon', 'RFvolc', 'RFsolar']:
    exec(VAR + '[:ind_attrib,0,:,:,...] = np.sum(' + VAR + '[:ind_attrib,:,:,:,...],1)')
    exec(VAR + '[:ind_attrib,1:,:,:,...] = 0')

# timeframed sectoral attribution
if (mod_sector == 'Time') & (300 < ind_final <= 310):
    for VAR in ['EFF', 'ECH4', 'EN2O'] + ['LUC', 'HARV', 'SHIFT'] + ['EHFC', 'EPFC',
                                                                     'EODS'] + ['ENOX',
                                                                                'ECO',
                                                                                'EVOC',
                                                                                'ESO2',
                                                                                'ENH3',
                                                                                'EOC',
                                                                                'EBC'] + [
                   'RFcon', 'RFvolc', 'RFsolar']:
        exec(VAR + '[150:200,:,:,1,:,...] = ' + VAR + '[150:200,:,:,0,:,...]')
        exec(VAR + '[200:210,:,:,2,:,...] = ' + VAR + '[200:210,:,:,0,:,...]')
        exec(VAR + '[210:220,:,:,3,:,...] = ' + VAR + '[210:220,:,:,0,:,...]')
        exec(VAR + '[220:230,:,:,4,:,...] = ' + VAR + '[220:230,:,:,0,:,...]')
        exec(VAR + '[230:240,:,:,5,:,...] = ' + VAR + '[230:240,:,:,0,:,...]')
        exec(VAR + '[240:250,:,:,6,:,...] = ' + VAR + '[240:250,:,:,0,:,...]')
        exec(VAR + '[250:260,:,:,7,:,...] = ' + VAR + '[250:260,:,:,0,:,...]')
        exec(VAR + '[260:265,:,:,8,:,...] = ' + VAR + '[260:265,:,:,0,:,...]')
        exec(VAR + '[265:270,:,:,9,:,...] = ' + VAR + '[265:270,:,:,0,:,...]')
        exec(VAR + '[270:275,:,:,10,:,...] = ' + VAR + '[270:275,:,:,0,:,...]')
        exec(VAR + '[275:280,:,:,11,:,...] = ' + VAR + '[275:280,:,:,0,:,...]')
        exec(VAR + '[280:285,:,:,12,:,...] = ' + VAR + '[280:285,:,:,0,:,...]')
        exec(VAR + '[285:290,:,:,13,:,...] = ' + VAR + '[285:290,:,:,0,:,...]')
        exec(VAR + '[290:292,:,:,14,:,...] = ' + VAR + '[290:292,:,:,0,:,...]')
        exec(VAR + '[292:294,:,:,15,:,...] = ' + VAR + '[292:294,:,:,0,:,...]')
        exec(VAR + '[294:296,:,:,16,:,...] = ' + VAR + '[294:296,:,:,0,:,...]')
        exec(VAR + '[296:298,:,:,17,:,...] = ' + VAR + '[296:298,:,:,0,:,...]')
        exec(VAR + '[298:300,:,:,18,:,...] = ' + VAR + '[298:300,:,:,0,:,...]')
        for t in range(300, ind_final + 1):
            exec(VAR + '[t,:,:,19+t-300,:,...] = ' + VAR + '[t,:,:,0,:,...]')
        exec(VAR + '[150:,:,:,0,:,...] = 0')
elif (mod_sector == 'TimeRCP') & (ind_final == 400):
    for VAR in ['EFF', 'ECH4', 'EN2O'] + ['LUC', 'HARV', 'SHIFT'] + ['EHFC', 'EPFC',
                                                                     'EODS'] + ['ENOX',
                                                                                'ECO',
                                                                                'EVOC',
                                                                                'ESO2',
                                                                                'ENH3',
                                                                                'EOC',
                                                                                'EBC'] + [
                   'RFcon', 'RFvolc', 'RFsolar']:
        exec(VAR + '[311:316,:,1,:,...] = ' + VAR + '[311:316,:,0,:,...]')
        exec(VAR + '[316:321,:,2,:,...] = ' + VAR + '[316:321,:,0,:,...]')
        exec(VAR + '[321:326,:,3,:,...] = ' + VAR + '[321:326,:,0,:,...]')
        exec(VAR + '[326:331,:,4,:,...] = ' + VAR + '[326:331,:,0,:,...]')
        exec(VAR + '[331:336,:,5,:,...] = ' + VAR + '[331:336,:,0,:,...]')
        exec(VAR + '[336:341,:,6,:,...] = ' + VAR + '[336:341,:,0,:,...]')
        exec(VAR + '[341:346,:,7,:,...] = ' + VAR + '[341:346,:,0,:,...]')
        exec(VAR + '[346:351,:,8,:,...] = ' + VAR + '[346:351,:,0,:,...]')
        exec(VAR + '[351:356,:,9,:,...] = ' + VAR + '[351:356,:,0,:,...]')
        exec(VAR + '[356:361,:,10,:,...] = ' + VAR + '[356:361,:,0,:,...]')
        exec(VAR + '[361:366,:,11,:,...] = ' + VAR + '[361:366,:,0,:,...]')
        exec(VAR + '[366:371,:,12,:,...] = ' + VAR + '[366:371,:,0,:,...]')
        exec(VAR + '[371:376,:,13,:,...] = ' + VAR + '[371:376,:,0,:,...]')
        exec(VAR + '[376:381,:,14,:,...] = ' + VAR + '[376:381,:,0,:,...]')
        exec(VAR + '[381:386,:,15,:,...] = ' + VAR + '[381:386,:,0,:,...]')
        exec(VAR + '[386:391,:,16,:,...] = ' + VAR + '[386:391,:,0,:,...]')
        exec(VAR + '[391:396,:,17,:,...] = ' + VAR + '[391:396,:,0,:,...]')
        exec(VAR + '[396:401,:,18,:,...] = ' + VAR + '[396:401,:,0,:,...]')
        exec(VAR + '[311:,:,0,:,...] = 0')











################################################################################
# OSCAR LOADP
################################################################################


print 'LOADING: PARAMETERS'

##################################################
#   1. CARBON DIOXIDE
##################################################

# ===============
# 1.1. ATMOSPHERE
# ===============

# conversion of CO2 from {ppm} to {GtC}
alpha_CO2 = 0.1765 * np.array([12.0], dtype=dty)

# historic CO2 from IPCC-AR5 {ppm}
# from [IPCC WG1, 2013] annexe 2
CO2_ipcc = np.ones([311 + 1], dtype=dty) * np.nan
TMP = np.array([line for line in csv.reader(
    open('data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1750-2011.CO2.csv', 'r'))],
               dtype=dty)
CO2_ipcc[50:] = TMP[:, 0]
CO2_0 = np.array([CO2_ipcc[50]], dtype=dty)

# historic CO2 from CMIP5 {ppm}
# from [Meinshausen et al., 2011]
CO2_cmip5 = np.ones([305 + 1], dtype=dty) * np.nan
TMP = np.array([line for line in csv.reader(
    open('data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005.CO2.csv', 'r'))], dtype=dty)
CO2_cmip5[65:] = TMP[:, 0]

# historic CO2 from NOAA {ppm}
# from the website
CO2_noaa_ml = np.ones([314 + 1], dtype=dty) * np.nan
TMP = np.array([line for line in csv.reader(
    open('data/HistAtmo_NOAA-ESRL/#DATA.HistAtmo_NOAA-ESRL.1959-2014.CO2_maunaloa.csv',
         'r'))], dtype=dty)
CO2_noaa_ml[259:] = TMP[:, 0]
CO2_noaa_gl = np.ones([314 + 1], dtype=dty) * np.nan
TMP = np.array([line for line in csv.reader(
    open('data/HistAtmo_NOAA-ESRL/#DATA.HistAtmo_NOAA-ESRL.1980-2014.CO2_global.csv',
         'r'))], dtype=dty)
CO2_noaa_gl[280:] = TMP[:, 0]

# historic CO2 from Law Dome ice cores {ppm}
# from [Etheridge et al., 1996] and [MacFarling Meure et al., 2006]
CO2_lawdome = np.array([line for line in csv.reader(
    open('data/HistAtmo_NOAA-NCDC/#DATA.HistAtmo_NOAA-NCDC.(IceCores).CO2_lawdome.csv',
         'r'))][1:], dtype=dty)

# load RCP concentrations {ppm}
# from [Meinshausen et al., 2011]
CO2_rcp = np.ones([800 + 1, 6], dtype=dty) * np.nan
n = -1
for rcp in ['rcp26', 'rcp45', 'rcp60', 'rcp85', 'rcp45to26', 'rcp60to45']:
    n += 1
    TMP = np.array([line for line in csv.reader(
        open('data/Scenario_ECP/#DATA.Scenario_ECP.2000-2500.' + rcp + '_CO2.csv', 'r'))],
                   dtype=dty)
    CO2_rcp[300:, n] = TMP[:, 0]

# global CO2 historic flux {GtC/yr}
# from [Le Quere et al., 2015]
TMP = np.array([line for line in csv.reader(
    open('data/Historic_GCP/#DATA.Historic_GCP.1959-2014_(5flx).budget.csv', 'r'))][1:],
               dtype=dty)
n = -1
for VAR in ['EFF', 'ELUC', 'd_CO2', 'OSNK', 'LSNK']:
    n += 1
    exec(VAR + '_gcp = np.ones([314+1], dtype=dty) * np.nan')
    exec(VAR + '_gcp[259:] = TMP[:,n]')
OSNK_gcp *= -1
LSNK_gcp *= -1

# ==========
# 1.2. OCEAN
# ==========

# ----------------
# 1.2.1. Structure
# ----------------

# parameters of the surface layer
# from HILDA and other models comparison [Joos et al., 1996]
# gas-exchange coefficient {/yr}, area {m2}, depth {m} and temperature {degC}
if (mod_OSNKstruct == 'HILDA'):
    v_fg = np.array([1 / 9.06], dtype=dty)
    A_ocean = np.array([3.62E14], dtype=dty)
    mld_0 = np.array([75], dtype=dty)
    sst_0 = np.array([18.2], dtype=dty)
elif (mod_OSNKstruct == 'BD-model'):
    v_fg = np.array([1 / 7.80], dtype=dty)
    A_ocean = np.array([3.62E14], dtype=dty)
    mld_0 = np.array([75], dtype=dty)
    sst_0 = np.array([17.7], dtype=dty)
elif (mod_OSNKstruct == '2D-model'):
    v_fg = np.array([1 / 7.46], dtype=dty)
    A_ocean = np.array([3.54E14], dtype=dty)
    mld_0 = np.array([50], dtype=dty)
    sst_0 = np.array([18.3], dtype=dty)
elif (mod_OSNKstruct == '3D-model'):
    v_fg = np.array([1 / 7.66], dtype=dty)
    A_ocean = np.array([3.55E14], dtype=dty)
    mld_0 = np.array([50.9], dtype=dty)
    sst_0 = np.array([17.7], dtype=dty)

# sizes and timescales of surface subdivisions for transportation to deep ocean
# adapted from the models comparison by [Joos et al., 1996]
nb_obox = 8
p_circ = np.zeros([nb_obox], dtype=dty)
tau_circ = np.ones([nb_obox], dtype=dty) * np.inf
if (mod_OSNKstruct == 'HILDA'):
    p_circ[1:1 + 6] = np.array([0.24278, 0.13963, 0.089318, 0.037820, 0.035549, 0.022936],
                               dtype=dty)
    p_circ[0] = 1 - np.sum(p_circ[1:])
    tau_circ[1:1 + 6] = np.array([1.2679, 5.2528, 18.601, 68.736, 232.30, np.inf],
                                 dtype=dty)
    tau_circ[0] = 1 / 3.
elif (mod_OSNKstruct == 'BD-model'):
    p_circ[1:1 + 7] = np.array(
        [0.16851, 0.11803, 0.076817, 0.050469, 0.010469, 0.031528, 0.019737], dtype=dty)
    p_circ[0] = 1 - np.sum(p_circ[:])
    tau_circ[1:1 + 7] = np.array([1.6388, 4.8702, 14.172, 43.506, 148.77, 215.71, np.inf],
                                 dtype=dty)
    tau_circ[0] = 1 / 3.
elif (mod_OSNKstruct == '2D-model'):
    p_circ[1:1 + 6] = np.array(
        [0.067380, 0.036608, 0.026994, 0.026933, 0.012456, 0.013691], dtype=dty)
    p_circ[0] = 1 - np.sum(p_circ[:])
    tau_circ[1:1 + 6] = np.array([10.515, 11.677, 38.946, 107.57, 331.54, np.inf],
                                 dtype=dty)
    tau_circ[0] = 1 / 2.
elif (mod_OSNKstruct == '3D-model'):
    p_circ[1:1 + 6] = np.array([0.70367, 0.24966, 0.066485, 0.038344, 0.019439, 0.014819],
                               dtype=dty)
    p_circ[1] = 1 - np.sum(p_circ[2:])
    tau_circ[1:1 + 6] = np.array([0.70177, 2.3488, 15.281, 65.359, 347.55, np.inf],
                                 dtype=dty)

# ----------------
# 1.2.2. Chemistry
# ----------------

# conversion factor for surface ocean{{mumol/kgsol}{m3/ppm}}
# after [Joos et al., 1996]
alpha_sol = np.array([1.722e17], dtype=dty)
alpha_dic = alpha_sol / alpha_CO2 / A_ocean / mld_0


# fonction for converting concentration to mass
def f_dic(D_CSURF, D_mld):
    D_dic = alpha_dic * D_CSURF / (1 + D_mld / mld_0)
    return np.array(D_dic, dtype=dty)


def df_dic(D_CSURF, D_mld):
    df_1 = alpha_dic / (1 + D_mld / mld_0) * D_CSURF
    df_2 = -alpha_dic * D_CSURF / mld_0 / (1 + D_mld / mld_0) ** 2 * D_mld
    df_tot = df_1 + df_2
    return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot)]


# emulation of the carbonate chemistry {ppm}
# after [Lewis and Wallas, 1998] and the CSIRO CTR047
if (mod_OSNKchem == 'CO2SysPade'):

    dic_0 = 30015.6 * (1 - 0.0226536 * (sst_0 - 15) + 0.000167105 * (sst_0 - 15) ** 2) * (
                CO2_0 / 380) / (1 + 13.4574 * (
                1 - 0.019829 * (sst_0 - 15) + 0.000113872 * (sst_0 - 15) ** 2) * (
                                            CO2_0 / 380) - 0.243121 * (1 + 0.000443511 * (
                sst_0 - 15) - 0.000473227 * (sst_0 - 15) ** 2) * (CO2_0 / 380) ** 2)
    CSURF_0 = p_circ * dic_0 / alpha_dic


    def f_pCO2(D_dic, D_sst):
        a0_0 = 30015.6 * (1 - 0.0226536 * (sst_0 - 15) + 0.000167105 * (sst_0 - 15) ** 2)
        a1_0 = 13.4574 * (1 - 0.019829 * (sst_0 - 15) + 0.000113872 * (sst_0 - 15) ** 2)
        a2_0 = -0.243121 * (
                    1 + 0.000443511 * (sst_0 - 15) - 0.000473227 * (sst_0 - 15) ** 2)
        a0 = 30015.6 * (1 - 0.0226536 * (sst_0 + D_sst - 15) + 0.000167105 * (
                    sst_0 + D_sst - 15) ** 2)
        a1 = 13.4574 * (1 - 0.019829 * (sst_0 + D_sst - 15) + 0.000113872 * (
                    sst_0 + D_sst - 15) ** 2)
        a2 = -0.243121 * (1 + 0.000443511 * (sst_0 + D_sst - 15) - 0.000473227 * (
                    sst_0 + D_sst - 15) ** 2)
        D_pCO2 = 380 * (a0 - a1 * (dic_0 + D_dic) - np.sqrt(
            (a0 - a1 * (dic_0 + D_dic)) ** 2 - 4 * a2 * (dic_0 + D_dic) ** 2)) / (
                             2 * a2 * (dic_0 + D_dic)) - 380 * (
                             a0_0 - a1_0 * dic_0 - np.sqrt(
                         (a0_0 - a1_0 * dic_0) ** 2 - 4 * a2_0 * dic_0 ** 2)) / (
                             2 * a2_0 * dic_0)
        return np.array(D_pCO2, dtype=dty)


    def df_pCO2_ddic(D_dic, D_sst):
        a0 = 30015.6 * (1 - 0.0226536 * (sst_0 + D_sst - 15) + 0.000167105 * (
                    sst_0 + D_sst - 15) ** 2)
        a1 = 13.4574 * (1 - 0.019829 * (sst_0 + D_sst - 15) + 0.000113872 * (
                    sst_0 + D_sst - 15) ** 2)
        a2 = -0.243121 * (1 + 0.000443511 * (sst_0 + D_sst - 15) - 0.000473227 * (
                    sst_0 + D_sst - 15) ** 2)
        rac = -np.sqrt((a0 - a1 * (dic_0 + D_dic)) ** 2 - 4 * a2 * (dic_0 + D_dic) ** 2)
        D_pCO2 = 380 * ((-2 * a0 * a1 * a2 * (dic_0 + D_dic) + (
                    2 * a2 * a1 ** 2 - 8 * a2 ** 2) * (
                                     dic_0 + D_dic) ** 2) / rac - 2 * a0 * a2 - 2 * a2 * rac) / (
                             2 * a2 * (dic_0 + D_dic)) ** 2
        return np.array(D_pCO2, dtype=dty)


    def df_pCO2_dsst(D_dic, D_sst):
        a0 = 30015.6 * (1 - 0.0226536 * (sst_0 + D_sst - 15) + 0.000167105 * (
                    sst_0 + D_sst - 15) ** 2)
        a1 = 13.4574 * (1 - 0.019829 * (sst_0 + D_sst - 15) + 0.000113872 * (
                    sst_0 + D_sst - 15) ** 2)
        a2 = -0.243121 * (1 + 0.000443511 * (sst_0 + D_sst - 15) - 0.000473227 * (
                    sst_0 + D_sst - 15) ** 2)
        da0 = 30015.6 * (-0.0226536 + 2 * 0.000167105 * (sst_0 + D_sst - 15))
        da1 = 13.4574 * (-0.019829 + 2 * 0.000113872 * (sst_0 + D_sst - 15))
        da2 = -0.243121 * (0.000443511 - 2 * 0.000473227 * (sst_0 + D_sst - 15))
        rac = -np.sqrt((a0 - a1 * (dic_0 + D_dic)) ** 2 - 4 * a2 * (dic_0 + D_dic) ** 2)
        D_pCO2 = 380 * (2 * (a2 * da0 - a0 * da2) * (dic_0 + D_dic) + 2 * (
                    a1 * da2 - da1 * a2) * (dic_0 + D_dic) ** 2 + (
                                    2 * a2 * da0 * a0 * (dic_0 + D_dic) - 2 * a2 * (
                                        a0 * da1 + a1 * da0) * (
                                                dic_0 + D_dic) ** 2 + 2 * a2 * (
                                                da1 * a1 - 2 * da2) * (
                                                dic_0 + D_dic) ** 3) / rac - 2 * da2 * (
                                    dic_0 + D_dic) * rac) / (
                             2 * a2 * (dic_0 + D_dic)) ** 2
        return np.array(D_pCO2, dtype=dty)


    def df_pCO2(D_dic, D_sst):
        df_1 = df_pCO2_ddic(D_dic, D_sst) * D_dic
        df_2 = df_pCO2_dsst(D_dic, D_sst) * D_sst
        df_tot = df_1 + df_2
        return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot)]

# after [Lewis and Wallas, 1998] and the CSIRO CTR047
elif (mod_OSNKchem == 'CO2SysPower'):

    dic_0 = 2160.156 * (
                1 - 0.00345063 * (sst_0 - 15) - 0.0000250016 * (sst_0 - 15) ** 2) * (
                        (CO2_0 / 380) - 0.318665 * (
                            1 - 0.00151292 * (sst_0 - 15) - 0.000198978 * (
                                sst_0 - 15) ** 2)) ** (0.0595961 * (
                1 + 0.0200328 * (sst_0 - 15) + 0.000192084 * (sst_0 - 15) ** 2))
    CSURF_0 = p_circ * dic_0 / alpha_dic


    def f_pCO2(D_dic, D_sst):
        p0_0 = 2160.156 * (
                    1 - 0.00345063 * (sst_0 - 15) - 0.0000250016 * (sst_0 - 15) ** 2)
        p1_0 = 0.0595961 * (
                    1 + 0.0200328 * (sst_0 - 15) + 0.000192084 * (sst_0 - 15) ** 2)
        p2_0 = 0.318665 * (
                    1 - 0.00151292 * (sst_0 - 15) - 0.000198978 * (sst_0 - 15) ** 2)
        p0 = 2160.156 * (1 - 0.00345063 * (sst_0 + D_sst - 15) - 0.0000250016 * (
                    sst_0 + D_sst - 15) ** 2)
        p1 = 0.0595961 * (1 + 0.0200328 * (sst_0 + D_sst - 15) + 0.000192084 * (
                    sst_0 + D_sst - 15) ** 2)
        p2 = 0.318665 * (1 - 0.00151292 * (sst_0 + D_sst - 15) - 0.000198978 * (
                    sst_0 + D_sst - 15) ** 2)
        D_pCO2 = 380 * (p2 + ((dic_0 + D_dic) / p0) ** (1 / p1)) - 380 * (
                    p2_0 + (dic_0 / p0_0) ** (1 / p1_0))
        return np.array(D_pCO2, dtype=dty)


    def df_pCO2_ddic(D_dic, D_sst):
        p0 = 2160.156 * (1 - 0.00345063 * (sst_0 + D_sst - 15) - 0.0000250016 * (
                    sst_0 + D_sst - 15) ** 2)
        p1 = 0.0595961 * (1 + 0.0200328 * (sst_0 + D_sst - 15) + 0.000192084 * (
                    sst_0 + D_sst - 15) ** 2)
        p2 = 0.318665 * (1 - 0.00151292 * (sst_0 + D_sst - 15) - 0.000198978 * (
                    sst_0 + D_sst - 15) ** 2)
        D_pCO2 = 380 * (1 / p0 / p1) * ((dic_0 + D_dic) / p0) ** (1 / p1 - 1)
        return np.array(D_pCO2, dtype=dty)


    def df_pCO2_dsst(D_dic, D_sst):
        p0 = 2160.156 * (1 - 0.00345063 * (sst_0 + D_sst - 15) - 0.0000250016 * (
                    sst_0 + D_sst - 15) ** 2)
        p1 = 0.0595961 * (1 + 0.0200328 * (sst_0 + D_sst - 15) + 0.000192084 * (
                    sst_0 + D_sst - 15) ** 2)
        p2 = 0.318665 * (1 - 0.00151292 * (sst_0 + D_sst - 15) - 0.000198978 * (
                    sst_0 + D_sst - 15) ** 2)
        dp0 = 2160.156 * (-0.00345063 - 2 * 0.0000250016 * (sst_0 + D_sst - 15))
        dp1 = 0.0595961 * (0.0200328 + 2 * 0.000192084 * (sst_0 + D_sst - 15))
        dp2 = 0.318665 * (-0.00151292 - 2 * 0.000198978 * (sst_0 + D_sst - 15))
        D_pCO2 = 380 * (dp2 + (
                    (-dp1 / p1 / p1) * np.log((dic_0 + D_dic) / p0) - (dp0 / p0 / p1)) * (
                                    (dic_0 + D_dic) / p0) ** (1 / p1))
        return np.array(D_pCO2, dtype=dty)


    def df_pCO2(D_dic, D_sst):
        df_1 = df_pCO2_ddic(D_dic, D_sst) * D_dic
        df_2 = df_pCO2_dsst(D_dic, D_sst) * D_sst
        df_tot = df_1 + df_2
        return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot)]

# ------------
# 1.2.3. CMIP5
# ------------

# period of CMIP5 simulations and depth of MLD per model
prd = {'ctrl': '451yr', 'hist': '1850-2005', 'rcp85': '2006-2300'}
lng = {'ctrl': 451, 'hist': 156, 'rcp85': 295}
mld = {'HILDA': '75m', 'BD-model': '75m', '2D-model': '50m', '3D-model': '50m9'}

# load pre-processed CMIP5 results for specified model
if (mod_OSNKtrans != ''):
    for sim in ['ctrl', 'hist', 'rcp85']:
        for VAR in ['FCO2', 'FEXP']:
            if os.path.isfile('data/Ocean_CMIP5/#DATA.Ocean_' + mod_OSNKtrans + '.' + prd[
                sim] + '_18x10lat.' + sim + '_' + VAR + '.csv'):
                TMP = np.array([line for line in csv.reader(open(
                    'data/Ocean_CMIP5/#DATA.Ocean_' + mod_OSNKtrans + '.' + prd[
                        sim] + '_18x10lat.' + sim + '_' + VAR + '.csv', 'r'))], dtype=dty)
                exec(VAR + '_' + sim + ' = np.sum(TMP,1)')
            else:
                exec(VAR + '_' + sim + ' = np.zeros([lng[sim]], dtype=dty)')
        for var in ['tos', 'dpCO2', 'sic', 'mld']:
            if os.path.isfile('data/Ocean_CMIP5/#DATA.Ocean_' + mod_OSNKtrans + '.' + prd[
                sim] + '_18x10lat.' + sim + '_' + var + '.csv'):
                TMP = np.array([line for line in csv.reader(open(
                    'data/Ocean_CMIP5/#DATA.Ocean_' + mod_OSNKtrans + '.' + prd[
                        sim] + '_18x10lat.' + sim + '_' + var + '.csv', 'r'))], dtype=dty)
                TMP2 = np.array([line for line in csv.reader(open(
                    'data/Ocean_CMIP5/#DATA.Ocean_' + mod_OSNKtrans + '.451yr_18x10lat.SURF.csv',
                    'r'))], dtype=dty)
                exec(
                    var + '_' + sim + ' = np.sum(TMP*TMP2[:len(TMP)],1)/np.sum(TMP2[:len(TMP)],1)')
            else:
                exec(var + '_' + sim + ' = np.zeros([lng[sim]], dtype=dty)')
        for var in ['amoc']:
            if os.path.isfile('data/Ocean_CMIP5/#DATA.Ocean_' + mod_OSNKtrans + '.' + prd[
                sim] + '.' + sim + '_' + var + '.csv'):
                TMP = np.array([line for line in csv.reader(open(
                    'data/Ocean_CMIP5/#DATA.Ocean_' + mod_OSNKtrans + '.' + prd[
                        sim] + '.' + sim + '_' + var + '.csv', 'r'))], dtype=dty)
                exec(var + '_' + sim + ' = TMP[:,0]')
            else:
                exec(var + '_' + sim + ' = np.zeros([lng[sim]], dtype=dty)')
    # aggregate all experiments
    if (mod_OSNKtrans != ''):
        for VAR in ['FCO2', 'FEXP'] + ['tos', 'dpCO2', 'sic', 'mld'] + ['amoc']:
            exec(VAR + '_all = np.array(list(' + VAR + '_hist)+list(' + VAR + '_rcp85))')

# definition of parameters
# for sea ice concentration reduction {pp}&{./K}
Alpha_sic = np.array([0], dtype=dty)
gamma_sic = np.array([0], dtype=dty)
# for sea surface exchange rate {.}&{./K}
gamma_fg = np.array([0], dtype=dty)
# for mixing layer stratification {.}&{./K}
alpha_mld = np.array([0], dtype=dty)
gamma_mld = np.array([0], dtype=dty)
# for overturning circulation slowdown {.}&{./K}
alpha_amoc = np.array([0], dtype=dty)
gamma_amoc = np.array([0], dtype=dty)
# for biological pump reduction decrease {GtC/yr}&{./K}
Alpha_exp = np.array([0], dtype=dty)
gamma_exp = np.array([0], dtype=dty)

# fit of parameters
# sic
if (mod_OSNKtrans != ''):
    obj = sic_all


    def err(var):
        carb = 1
        clim = np.mean(sic_ctrl) * np.exp(
            1 - np.exp(var[0] * (tos_all - np.mean(tos_ctrl))))
        return np.sum((obj - carb * clim) ** 2)


    [gamma_sic[0]] = fmin(err, [0], disp=False)
    Alpha_sic[0] = np.mean(sic_ctrl)
# fg
if (mod_OSNKtrans != ''):
    diff = FCO2_all - np.mean(FCO2_ctrl)


    def err(var):
        carb = 1
        clim = var[0] * (1 - sic_all) * (dpCO2_all - np.mean(dpCO2_ctrl)) * (
                    1 + var[1] * (tos_all - np.mean(tos_ctrl)))
        return np.sum((diff - carb * clim) ** 2)


    [TMP, gamma_fg[0]] = fmin(err, [-1, 0], disp=False)
# mld
if (mod_OSNKtrans != ''):
    ratio = mld_all / np.mean(mld_ctrl)


    def err(var):
        carb = 1
        clim = (1 - var[0]) + var[0] * np.exp(var[1] * (tos_all - np.mean(tos_ctrl)))
        return np.sum((ratio - carb * clim) ** 2)


    [alpha_mld[0], gamma_mld[0]] = fmin(err, [0.5, 0], disp=False)
# amoc
if (mod_OSNKtrans != ''):
    ratio = amoc_all / np.mean(amoc_ctrl)


    def err(var):
        carb = 1
        clim = (1 - var[0]) + var[0] * np.exp(
            1 - np.exp(var[1] * (tos_all - np.mean(tos_ctrl))))
        return np.sum((ratio - carb * clim) ** 2)


    [alpha_amoc[0], gamma_amoc[0]] = fmin(err, [0.5, 0], disp=False)
# exp
if (mod_OSNKtrans != ''):
    diff = FEXP_all - np.mean(FEXP_ctrl)


    def err(var):
        carb = 1
        clim = var[0] * (1 - np.exp(var[1] * (tos_all - np.mean(tos_ctrl))))
        return np.sum((diff - carb * clim) ** 2)


    [Alpha_exp[0], gamma_exp[0]] = fmin(err, [-1, 0], disp=False)

# =========
# 1.3. LAND
# =========

# ----------------
# 1.3.1. Functions
# ----------------

# reference beta-factor of NPP {.}
# from FACE experiments [Norby et al., 2005]
beta_ref = np.array([0.60], dtype=dty)

# reference Q10 of HR {.}
# from [Mahecha et al., 2010]
Q10_ref = np.array([1.4], dtype=dty)

# expression of NPP function
# logarithmic formulation [Friedlingstein et al., 1995]
if (mod_LSNKnpp == 'log'):
    def f_npp(D_CO2, D_lst, D_lyp):
        D_npp = (1 + beta_npp * np.log(1 + D_CO2 / CO2_0)) * (
                    1 + gamma_nppT * D_lst + gamma_nppP * D_lyp) - 1
        return np.array(D_npp, dtype=dty)


    def df_npp_dCO2(D_CO2, D_lst, D_lyp):
        D_npp = beta_npp * (1 + gamma_nppT * D_lst + gamma_nppP * D_lyp) / (CO2_0 + D_CO2)
        return np.array(D_npp, dtype=dty)


    def df_npp_dlst(D_CO2, D_lst, D_lyp):
        D_npp = gamma_nppT * (1 + beta_npp * np.log(1 + D_CO2 / CO2_0))
        return np.array(D_npp, dtype=dty)


    def df_npp_dlyp(D_CO2, D_lst, D_lyp):
        D_npp = gamma_nppP * (1 + beta_npp * np.log(1 + D_CO2 / CO2_0))
        return np.array(D_npp, dtype=dty)


    def df_npp(D_CO2, D_lst, D_lyp):
        df_1 = df_npp_dCO2(D_CO2, D_lst, D_lyp) * D_CO2
        df_2 = df_npp_dlst(D_CO2, D_lst, D_lyp) * D_lst
        df_3 = df_npp_dlyp(D_CO2, D_lst, D_lyp) * D_lyp
        df_tot = df_1 + df_2 + df_3
        return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot),
                np.nan_to_num(df_3 / df_tot)]

# hyperbolic formulation [Friedlingstein et al., 1995]
elif (mod_LSNKnpp == 'hyp'):
    def f_npp(D_CO2, D_lst, D_lyp):
        K2 = ((2 * CO2_0 - CO2_comp) - beta_npp0 * (CO2_0 - CO2_comp)) / (
                    (beta_npp0 - 1) * (CO2_0 - CO2_comp) * (2 * CO2_0 - CO2_comp))
        K1 = (1 + K2 * (CO2_0 - CO2_comp)) / (CO2_0 - CO2_comp)
        D_npp = K1 * (CO2_0 + D_CO2 - CO2_comp) / (
                    1 + K2 * (CO2_0 + D_CO2 - CO2_comp)) * (
                            1 + gamma_nppT * D_lst + gamma_nppP * D_lyp) - 1
        return np.array(D_npp, dtype=dty)


    def df_npp_dCO2(D_CO2, D_lst, D_lyp):
        K2 = ((2 * CO2_0 - CO2_comp) - beta_npp0 * (CO2_0 - CO2_comp)) / (
                    (beta_npp0 - 1) * (CO2_0 - CO2_comp) * (2 * CO2_0 - CO2_comp))
        K1 = (1 + K2 * (CO2_0 - CO2_comp)) / (CO2_0 - CO2_comp)
        D_npp = (1 + gamma_nppT * D_lst + gamma_nppP * D_lyp) * K1 / (
                    1 + K2 * (CO2_0 + D_CO2 - CO2_comp)) ** 2
        return np.array(D_npp, dtype=dty)


    def df_npp_dlst(D_CO2, D_lst, D_lyp):
        K2 = ((2 * CO2_0 - CO2_comp) - beta_npp0 * (CO2_0 - CO2_comp)) / (
                    (beta_npp0 - 1) * (CO2_0 - CO2_comp) * (2 * CO2_0 - CO2_comp))
        K1 = (1 + K2 * (CO2_0 - CO2_comp)) / (CO2_0 - CO2_comp)
        D_npp = gamma_nppT * K1 * (CO2_0 + D_CO2 - CO2_comp) / (
                    1 + K2 * (CO2_0 + D_CO2 - CO2_comp))
        return np.array(D_npp, dtype=dty)


    def df_npp_dlyp(D_CO2, D_lst, D_lyp):
        K2 = ((2 * CO2_0 - CO2_comp) - beta_npp0 * (CO2_0 - CO2_comp)) / (
                    (beta_npp0 - 1) * (CO2_0 - CO2_comp) * (2 * CO2_0 - CO2_comp))
        K1 = (1 + K2 * (CO2_0 - CO2_comp)) / (CO2_0 - CO2_comp)
        D_npp = gamma_nppP * K1 * (CO2_0 + D_CO2 - CO2_comp) / (
                    1 + K2 * (CO2_0 + D_CO2 - CO2_comp))
        return np.array(D_npp, dtype=dty)


    def df_npp(D_CO2, D_lst, D_lyp):
        df_1 = df_npp_dCO2(D_CO2, D_lst, D_lyp) * D_CO2
        df_2 = df_npp_dlst(D_CO2, D_lst, D_lyp) * D_lst
        df_3 = df_npp_dlyp(D_CO2, D_lst, D_lyp) * D_lyp
        df_tot = df_1 + df_2 + df_3
        return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot),
                np.nan_to_num(df_3 / df_tot)]

# expression of RH function
# exponential formulation [Tuomi et al., 2008]
if (mod_LSNKrho == 'exp'):
    def f_rho(D_lst, D_lyp):
        D_rho = np.exp(gamma_rhoT * D_lst + gamma_rhoP * D_lyp) - 1
        return np.array(D_rho, dtype=dty)


    def df_rho_dlst(D_lst, D_lyp):
        D_rho = gamma_rhoT * np.exp(gamma_rhoT * D_lst + gamma_rhoP * D_lyp)
        return np.array(D_rho, dtype=dty)


    def df_rho_dlyp(D_lst, D_lyp):
        D_rho = gamma_rhoP * np.exp(gamma_rhoT * D_lst + gamma_rhoP * D_lyp)
        return np.array(D_rho, dtype=dty)


    def df_rho(D_lst, D_lyp):
        df_1 = df_rho_dlst(D_lst, D_lyp) * D_lst
        df_2 = df_rho_dlyp(D_lst, D_lyp) * D_lyp
        df_tot = df_1 + df_2
        return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot)]

# gaussian formulation [Tuomi et al., 2008]
elif (mod_LSNKrho == 'gauss'):
    def f_rho(D_lst, D_lyp):
        D_rho = np.exp(
            gamma_rhoT1 * D_lst + gamma_rhoT2 * D_lst ** 2 + gamma_rhoP * D_lyp) - 1
        return np.array(D_rho, dtype=dty)


    def df_rho_dlst(D_lst, D_lyp):
        D_rho = (gamma_rhoT1 + 2 * gamma_rhoT2 * D_lst) * np.exp(
            gamma_rhoT1 * D_lst + gamma_rhoT2 * D_lst ** 2 + gamma_rhoP * D_lyp)
        return np.array(D_rho, dtype=dty)


    def df_rho_dlyp(D_lst, D_lyp):
        D_rho = gamma_rhoP * np.exp(
            gamma_rhoT1 * D_lst + gamma_rhoT2 * D_lst ** 2 + gamma_rhoP * D_lyp)
        return np.array(D_rho, dtype=dty)


    def df_rho(D_lst, D_lyp):
        df_1 = df_rho_dlst(D_lst, D_lyp) * D_lst
        df_2 = df_rho_dlyp(D_lst, D_lyp) * D_lyp
        df_tot = df_1 + df_2
        return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot)]

# -------------
# 1.3.2. TRENDY
# -------------

# basic biomes of aggregation
bio = ['des', 'for', 'shr', 'gra', 'cro', 'pas', 'urb']

# load pre-processed TRENDY results for specified model
# data related to preindustrial C-cycle
for VAR in ['AREA', 'CSOIL', 'CLITTER', 'CVEG', 'NPP', 'RH']:
    exec(VAR + '_pi = np.zeros([nb_regionI,len(bio)], dtype=dty)')
    if os.path.isfile(
            'data/Land_TRENDYv2/#DATA.Land_' + mod_LSNKpreind + '.1910s_114reg1_7bio.' + VAR + '.csv'):
        TMP = np.array([line for line in csv.reader(open(
            'data/Land_TRENDYv2/#DATA.Land_' + mod_LSNKpreind + '.1910s_114reg1_7bio.' + VAR + '.csv',
            'r'))], dtype=dty)
        for i in range(1, 114 + 1):
            exec(VAR + '_pi[regionI_index[i],:] += TMP[i-1,:]')

# complete with arbitrary values
# if no litter given by model
if (np.sum(CLITTER_pi) == 0):
    CLITTER_pi = 0.05 * CSOIL_pi
    CSOIL_pi = 0.95 * CSOIL_pi
# if no shrubs in the model
if (nb_biome > 1) & (np.sum(AREA_pi[:, bio.index("shr")]) == 0) & (mod_biomeSHR == 'SHR'):
    AREA_pi[:, bio.index("shr")] = AREA_pi[:, bio.index("gra")]
    for VAR in ['CSOIL', 'CLITTER', 'CVEG', 'NPP', 'RH']:
        exec(
            VAR + '_pi[:,bio.index("shr")] = 0.85*' + VAR + '_pi[:,bio.index("gra")] + 0.15*' + VAR + '_pi[:,bio.index("for")]*AREA_pi[:,bio.index("gra")]/AREA_pi[:,bio.index("for")]')
    TMP = (AREA_pi[:, :, bio.index("for")] == 0)
    exec(
        VAR + '_pi[:,:,bio.index("shr")][TMP] = ' + VAR + '_pi[:,:,bio.index("gra")][TMP]')
# if no crops in the model
if (nb_biome > 1) & (np.sum(AREA_pi[:, bio.index("cro")]) == 0):
    for VAR in ['AREA', 'CSOIL', 'CLITTER', 'CVEG', 'NPP', 'RH']:
        exec(VAR + '_pi[:,bio.index("cro")] = ' + VAR + '_pi[:,bio.index("gra")]')
# if no pastures in the model
if (nb_biome > 1) & (np.sum(AREA_pi[:, bio.index("pas")]) == 0):
    AREA_pi[:, bio.index("pas")] = AREA_pi[:, bio.index("gra")]
    for VAR in ['CSOIL', 'CLITTER', 'CVEG', 'NPP', 'RH']:
        exec(
            VAR + '_pi[:,bio.index("pas")] = 0.60*' + VAR + '_pi[:,bio.index("gra")] + 0.40*' + VAR + '_pi[:,bio.index("des")]*AREA_pi[:,bio.index("gra")]/AREA_pi[:,bio.index("des")]')
        TMP = (AREA_pi[:, bio.index("des")] == 0)
        exec(
            VAR + '_pi[:,bio.index("pas")][TMP] = ' + VAR + '_pi[:,bio.index("gra")][TMP]')
# if no urban in the model
if (nb_biome > 1) & (np.sum(AREA_pi[:, bio.index("urb")]) == 0) & (
        (mod_biomeURB == 'URB') | mod_biomeV3):
    for VAR in ['AREA', 'CSOIL', 'CLITTER', 'CVEG', 'NPP', 'RH']:
        exec(VAR + '_pi[:,bio.index("urb")] = ' + VAR + '_pi[:,bio.index("des")]')

# data related to preindustrial fires
for VAR in ['EFIRE', 'CBURNT']:
    exec(VAR + '_pi = np.zeros([nb_regionI,len(bio)], dtype=dty)')
    if (mod_EFIREpreind != ''):
        if os.path.isfile(
                'data/BioBurn_TRENDYv2/#DATA.BioBurn_' + mod_EFIREpreind + '.1910s_114reg1_7bio.' + VAR + '.csv'):
            TMP = np.array([line for line in csv.reader(open(
                'data/BioBurn_TRENDYv2/#DATA.BioBurn_' + mod_EFIREpreind + '.1910s_114reg1_7bio.' + VAR + '.csv',
                'r'))], dtype=dty)
        for i in range(1, 114 + 1):
            exec(VAR + '_pi[regionI_index[i],:] += TMP[i-1,:]')

# complete with arbitrary values
# if no shrubs in the model
if (nb_biome > 1) & (np.sum(CBURNT_pi[:, bio.index("shr")]) == 0) & (
        mod_biomeSHR == 'SHR'):
    CBURNT_pi[:, bio.index("shr")] = CBURNT_pi[:, bio.index("gra")]
    for VAR in ['EFIRE']:
        exec(
            VAR + '_pi[:,bio.index("shr")] = 0.85*' + VAR + '_pi[:,bio.index("gra")] + 0.15*' + VAR + '_pi[:,bio.index("for")]*CBURNT_pi[:,bio.index("gra")]/CBURNT_pi[:,bio.index("for")]')
        TMP = (CBURNT_pi[:, :, bio.index("for")] == 0)
        exec(
            VAR + '_pi[:,:,bio.index("shr")][TMP] = ' + VAR + '_pi[:,:,bio.index("gra")][TMP]')
# if no crops in the model
if (nb_biome > 1) & (np.sum(CBURNT_pi[:, bio.index("cro")]) == 0):
    for VAR in ['EFIRE', 'CBURNT']:
        exec(VAR + '_pi[:,bio.index("cro")] = ' + VAR + '_pi[:,bio.index("gra")]')
# if no pastures in the model
if (nb_biome > 1) & (np.sum(CBURNT_pi[:, bio.index("pas")]) == 0):
    CBURNT_pi[:, bio.index("pas")] = CBURNT_pi[:, bio.index("gra")]
    for VAR in ['EFIRE']:
        exec(
            VAR + '_pi[:,bio.index("pas")] = 0.60*' + VAR + '_pi[:,bio.index("gra")] + 0.40*' + VAR + '_pi[:,bio.index("des")]*CBURNT_pi[:,bio.index("gra")]/CBURNT_pi[:,bio.index("des")]')
        TMP = (CBURNT_pi[:, bio.index("des")] == 0)
        exec(
            VAR + '_pi[:,bio.index("pas")][TMP] = ' + VAR + '_pi[:,bio.index("gra")][TMP]')
# if no urban in the model
if (nb_biome > 1) & (np.sum(CBURNT_pi[:, bio.index("urb")]) == 0) & (
        (mod_biomeURB == 'URB') | mod_biomeV3):
    for VAR in ['EFIRE', 'CBURNT']:
        exec(VAR + '_pi[:,bio.index("urb")] = ' + VAR + '_pi[:,bio.index("des")]')

# aggregate biomes
for VAR in ['AREA', 'CSOIL', 'CLITTER', 'CVEG', 'NPP', 'RH'] + ['EFIRE', 'CBURNT']:
    exec('TMP = ' + VAR + '_pi.copy()')
    exec(VAR + '_pi = np.zeros([nb_regionI,nb_biome], dtype=dty)')
    for b in range(len(bio)):
        exec(VAR + '_pi[:,biome_index[bio[b]]] += TMP[:,b]')

# ratio for metabolization of litter to soil {.}
# from [Foley, 1995]
k_met = np.array([0.3 / 0.7], dtype=dty)

# calculate preindustrial flux parameters
# net primary productivity {GtC/Mha/yr}
npp_0 = NPP_pi / AREA_pi
# rate of biomass mortality {./yr}
mu_0 = RH_pi / CVEG_pi
# rates of heterotrophic respiration {./yr}
rho_0 = RH_pi / (CLITTER_pi + CSOIL_pi)
rho1_0 = 1 / (1 + k_met) * RH_pi / CLITTER_pi
rho2_0 = k_met / (1 + k_met) * RH_pi / CSOIL_pi
# rate of fire ignition {./yr}
igni_0 = EFIRE_pi / CBURNT_pi
# [NaN]
for var in ['npp', 'mu', 'rho', 'rho1', 'rho2', 'igni']:
    exec(var + '_0[np.isnan(' + var + '_0)|np.isinf(' + var + '_0)] = 0')

# arbitrary adjustment of parameters
# lower preindustrial npp
npp_0 *= CO2_0 / np.mean(CO2_cmip5[201:231])
# crops: yearly harvest of 80% & protection against wildfires of 100%
if (nb_biome > 1):
    npp_0[:, biome_index["cro"]] *= 0.2
    mu_0[:, biome_index["cro"]] = 1.
    igni_0[:, biome_index["cro"]] *= 0.0
# pasts: harvest of 0% & protection against wildfires of 0%
if (nb_biome > 1):
    npp_0[:, biome_index["pas"]] *= 1.0
    igni_0[:, biome_index["pas"]] *= 1.0
# urban: protection against wildfires of 100%
if (nb_biome > 1) & ((mod_biomeURB == 'URB') | mod_biomeV3):
    igni_0[:, biome_index["urb"]] = 0.

# calculate preindustrial stocks {GtC/Mha}
cveg_0 = npp_0 / (mu_0 + igni_0)
csoil_0 = mu_0 / rho_0 * cveg_0
csoil1_0 = 1 / (1 + k_met) * mu_0 / rho1_0 * cveg_0
csoil2_0 = k_met * (rho1_0 / rho2_0) * csoil1_0
# [NaN]
for var in ['cveg', 'csoil', 'csoil1', 'csoil2']:
    exec(var + '_0[np.isnan(' + var + '_0)|np.isinf(' + var + '_0)] = 0')

# -----------------
# 1.3.3. Land-Cover
# -----------------

# load preindustrial land-cover {Mha}
AREA_0 = np.zeros([nb_regionI, nb_biome], dtype=dty)
TMP = np.array([line for line in csv.reader(open(
    'data/LandCover_' + data_LULCC + '/#DATA.LandCover_' + data_LULCC + '_' + mod_LSNKcover + '.1700_114reg1_7bio.AREA.csv',
    'r'))], dtype=dty)
for i in range(1, 114 + 1):
    for b in range(len(bio)):
        AREA_0[regionI_index[i], biome_index[bio[b]]] += TMP[i - 1, b]

# force true 1750 preindustrial (land-cover)
if PI_1750:
    LUC_tmp = np.zeros([50 + 1, nb_regionI, nb_biome, nb_biome], dtype=dty)
    bio = ['des', 'for', 'shr', 'gra', 'cro', 'pas', 'urb']
    # load LUH data
    if (data_LULCC[:3] == 'LUH'):
        for b1 in range(len(bio)):
            for b2 in range(len(bio)):
                if os.path.isfile(
                        'data/LandUse_' + data_LULCC + '/#DATA.LandUse_' + data_LULCC + '_' + mod_LSNKcover + '.1501-2015_114reg1.LUC_' +
                        bio[b1] + '2' + bio[b2] + '.csv'):
                    TMP = np.array([line for line in csv.reader(open(
                        'data/LandUse_LUH1/#DATA.LandUse_' + data_LULCC + '_' + mod_LSNKcover + '.1501-2015_114reg1.LUC_' +
                        bio[b1] + '2' + bio[b2] + '.csv', 'r'))], dtype=dty)
                    for i in range(1, 114 + 1):
                        LUC_tmp[1:50 + 1, regionI_index[i], biome_index[bio[b1]],
                        biome_index[bio[b2]]] += TMP[200:200 + 50, i - 1]
    # correct land-cover
    AREA_0 += np.sum(np.sum(LUC_tmp, 2), 0) - np.sum(np.sum(LUC_tmp, 3), 0)
    del LUC_tmp

# ------------
# 1.3.4. CMIP5
# ------------

# basic biomes of aggregation
bio = ['des', 'for', 'shr', 'gra', 'cro', 'pas', 'urb']

# load pre-processed CMIP5 results for specified model
# data related to C-cycle
for sim in ['ctrl', 'upct', 'fxcl', 'fdbk']:
    for VAR in ['AREA', 'CSOIL0', 'NPP', 'RH', 'FINPUT']:
        exec(VAR + '_' + sim + ' = np.zeros([140,nb_regionI,len(bio)], dtype=dty)')
        for b in range(len(bio)):
            if os.path.isfile(
                    'data/Land_CMIP5/#DATA.Land_' + mod_LSNKtrans + '.140yr_114reg1.' + sim + '_' + VAR + '_' +
                    bio[b] + '.csv'):
                TMP = np.array([line for line in csv.reader(open(
                    'data/Land_CMIP5/#DATA.Land_' + mod_LSNKtrans + '.140yr_114reg1.' + sim + '_' + VAR + '_' +
                    bio[b] + '.csv', 'r'))], dtype=dty)
                for i in range(1, 114 + 1):
                    exec(VAR + '_' + sim + '[:,regionI_index[i],b] += TMP[:,i-1]')
    for var in ['tas', 'pr']:
        exec(var + '_' + sim + ' = np.zeros([140,nb_regionI], dtype=dty)')
        TMP = np.array([line for line in csv.reader(open(
            'data/Land_CMIP5/#DATA.Land_' + mod_LSNKtrans + '.140yr_114reg1.' + sim + '_' + var + '.csv',
            'r'))], dtype=dty)
        TMP2 = np.array([line for line in csv.reader(open(
            'data/Land_CMIP5/#DATA.Land_' + mod_LSNKtrans + '.140yr_114reg1.AREA.csv',
            'r'))], dtype=dty)
        for i in range(1, 114 + 1):
            exec(var + '_' + sim + '[:,regionI_index[i]] += TMP[:,i-1]*TMP2[:,i-1]')
        TMP = np.zeros([140, nb_regionI], dtype=dty)
        for i in range(1, 114 + 1):
            exec('TMP[:,regionI_index[i]] += TMP2[:,i-1]')
        exec(var + '_' + sim + ' /= TMP')
        exec(
            var + '_' + sim + '[np.isnan(' + var + '_' + sim + ')|np.isinf(' + var + '_' + sim + ')] = 0')

# complete with arbitrary values
# if no shrubs in the model
if (nb_biome > 1) & (np.sum(AREA_ctrl[:, :, bio.index("shr")]) == 0) & (
        mod_biomeSHR == 'SHR'):
    for sim in ['ctrl', 'upct', 'fxcl', 'fdbk']:
        exec(
            'AREA_' + sim + '[:,:,bio.index("shr")] = AREA_' + sim + '[:,:,bio.index("gra")]')
        for VAR in ['CSOIL0', 'NPP', 'RH', 'FINPUT']:
            exec(
                VAR + '_' + sim + '[:,:,bio.index("shr")] = 0.85*' + VAR + '_' + sim + '[:,:,bio.index("gra")] + 0.15*' + VAR + '_' + sim + '[:,:,bio.index("for")]*AREA_' + sim + '[:,:,bio.index("gra")]/AREA_' + sim + '[:,:,bio.index("for")]')
            exec('TMP = (AREA_' + sim + '[:,:,bio.index("for")] == 0)')
            exec(
                VAR + '_' + sim + '[:,:,bio.index("shr")][TMP] = ' + VAR + '_' + sim + '[:,:,bio.index("gra")][TMP]')
# if no crops in the model
if (nb_biome > 1) & (np.sum(AREA_ctrl[:, :, bio.index("cro")]) == 0):
    for sim in ['ctrl', 'upct', 'fxcl', 'fdbk']:
        for VAR in ['AREA', 'CSOIL0', 'NPP', 'RH', 'FINPUT']:
            exec(
                VAR + '_' + sim + '[:,:,bio.index("cro")] = ' + VAR + '_' + sim + '[:,:,bio.index("gra")]')
# if no pastures in the model
if (nb_biome > 1) & (np.sum(AREA_ctrl[:, :, bio.index("pas")]) == 0):
    for sim in ['ctrl', 'upct', 'fxcl', 'fdbk']:
        exec(
            'AREA_' + sim + '[:,:,bio.index("pas")] = AREA_' + sim + '[:,:,bio.index("gra")]')
        for VAR in ['CSOIL0', 'NPP', 'RH', 'FINPUT']:
            exec(
                VAR + '_' + sim + '[:,:,bio.index("pas")] = 0.60*' + VAR + '_' + sim + '[:,:,bio.index("gra")] + 0.40*' + VAR + '_' + sim + '[:,:,bio.index("des")]*AREA_' + sim + '[:,:,bio.index("gra")]/AREA_' + sim + '[:,:,bio.index("des")]')
            exec('TMP = (AREA_' + sim + '[:,:,bio.index("des")] == 0)')
            exec(
                VAR + '_' + sim + '[:,:,bio.index("pas")][TMP] = ' + VAR + '_' + sim + '[:,:,bio.index("gra")][TMP]')
# if no urban in the model
if (nb_biome > 1) & (np.sum(AREA_ctrl[:, :, bio.index("urb")]) == 0) & (
        (mod_biomeURB == 'URB') | mod_biomeV3):
    for sim in ['ctrl', 'upct', 'fxcl', 'fdbk']:
        for VAR in ['AREA', 'CSOIL0', 'NPP', 'RH', 'FINPUT']:
            exec(
                VAR + '_' + sim + '[:,:,bio.index("urb")] = ' + VAR + '_' + sim + '[:,:,bio.index("des")]')

# aggregate biomes
for sim in ['ctrl', 'upct', 'fxcl', 'fdbk']:
    for VAR in ['AREA', 'CSOIL0', 'NPP', 'RH', 'FINPUT']:
        exec('TMP = ' + VAR + '_' + sim + '.copy()')
        exec(VAR + '_' + sim + ' = np.zeros([140,nb_regionI,nb_biome], dtype=dty)')
        for b in range(len(bio)):
            exec(VAR + '_' + sim + '[:,:,biome_index[bio[b]]] += TMP[:,:,b]')
# aggregate all experiments
CO2_all = np.array(list(CO2_cmip5[150] * 1.01 ** np.arange(140)) + list(
    CO2_cmip5[150] * 1.01 ** np.arange(140)) + list(CO2_cmip5[150] * np.ones(140)),
                   dtype=dty)
for VAR in ['AREA', 'CSOIL0', 'NPP', 'RH', 'FINPUT'] + ['tas', 'pr']:
    exec(
        VAR + '_all = np.array(list(' + VAR + '_upct)+list(' + VAR + '_fxcl)+list(' + VAR + '_fdbk), dtype=dty)')
# decadal means
for VAR in ['AREA', 'CSOIL0', 'NPP', 'RH', 'FINPUT'] + ['tas', 'pr'] + ['CO2']:
    exec(
        VAR + '_dec= np.zeros([3*(140-10)]+list(np.shape(' + VAR + '_all)[1:]), dtype=dty)')
    for t in range(130):
        exec(VAR + '_dec[t,...] = np.mean(' + VAR + '_all[t:t+10,...],0)')
        exec(VAR + '_dec[t+140-10,...] = np.mean(' + VAR + '_all[140+t:140+t+10,...],0)')
        exec(
            VAR + '_dec[t+2*140-2*10,...] = np.mean(' + VAR + '_all[2*140+t:2*140+t+10,...],0)')

# definition of parameters
# sensitivity of NPP to CO2 {.} or {.}&{ppm}
if (mod_LSNKnpp == 'log'):
    beta_npp = np.zeros([nb_regionI, nb_biome], dtype=dty)
elif (mod_LSNKnpp == 'hyp'):
    beta_npp0 = np.ones([nb_regionI, nb_biome], dtype=dty) * 1.000001
    CO2_comp = np.ones([nb_regionI, nb_biome], dtype=dty) * 1E18
# sensitivity of NPP to climate {/K}&{/mm}
gamma_nppT = np.zeros([nb_regionI, nb_biome], dtype=dty)
gamma_nppP = np.zeros([nb_regionI, nb_biome], dtype=dty)
# sensitivity of RH rate to soil input {.} (not used)
prim_rho = np.zeros([nb_regionI, nb_biome], dtype=dty)
# sensitivity of RH rate to climate {/K}&{/mm} or {/K}&{/K**2}&{/mm}
if (mod_LSNKrho == 'exp'):
    gamma_rhoT = np.zeros([nb_regionI, nb_biome], dtype=dty)
    gamma_rhoP = np.zeros([nb_regionI, nb_biome], dtype=dty)
elif (mod_LSNKrho == 'gauss'):
    gamma_rhoT1 = np.zeros([nb_regionI, nb_biome], dtype=dty)
    gamma_rhoT2 = np.zeros([nb_regionI, nb_biome], dtype=dty)
    gamma_rhoP = np.zeros([nb_regionI, nb_biome], dtype=dty)

# fit of parameters
for i in range(1, nb_regionI):
    for b in range(nb_biome):

        # for NPP
        ratio = (NPP_all / AREA_all)[:, i, b] / np.mean((NPP_ctrl / AREA_ctrl)[:, i, b])
        ratio_dec = (NPP_dec / AREA_dec)[:, i, b] / np.mean(
            (NPP_ctrl / AREA_ctrl)[:, i, b])
        if (mod_LSNKnpp == 'log'):
            def err(var):
                carb = 1 + np.abs(var[0]) * np.log(CO2_dec[:] / CO2_cmip5[150])
                clim = 1 + var[1] * (tas_dec[:, i] - np.mean(tas_ctrl[:, i]))
                return np.sum((ratio_dec - carb * clim) ** 2)


            [beta_npp[i, b], gamma_nppT[i, b]] = fmin(err, [beta_ref[0], 0], disp=False)
            beta_npp[i, b] = np.abs(beta_npp[i, b])


            def err(var):
                carb = 1 + beta_npp[i, b] * np.log(CO2_all[:] / CO2_cmip5[150])
                clim = 1 + gamma_nppT[i, b] * (tas_all[:, i] - np.mean(tas_ctrl[:, i])) + \
                       var[0] * (pr_all[:, i] - np.mean(pr_ctrl[:, i]))
                return np.sum((ratio - carb * clim) ** 2)


            [gamma_nppP[i, b]] = fmin(err, [0], disp=False)
        elif (mod_LSNKnpp == 'hyp'):
            def err(var):
                K2 = ((2 * CO2_cmip5[150] - np.abs(var[1])) - (1 + np.abs(var[0])) * (
                            CO2_cmip5[150] - np.abs(var[1]))) / (
                                 (1 + np.abs(var[0]) - 1) * (
                                     CO2_cmip5[150] - np.abs(var[1])) * (
                                             2 * CO2_cmip5[150] - np.abs(var[1])))
                K1 = (1 + K2 * (CO2_cmip5[150] - np.abs(var[1]))) / (
                            CO2_cmip5[150] - np.abs(var[1]))
                carb = (K1 * (CO2_dec[:] - np.abs(var[1]))) / (
                            1 + K2 * (CO2_dec[:] - np.abs(var[1])))
                clim = 1 + var[2] * (tas_dec[:, i] - np.mean(tas_ctrl[:, i]))
                return np.sum((ratio_dec - carb * clim) ** 2)


            [beta_npp0[i, b], CO2_comp[i, b], gamma_nppT[i, b]] = fmin(err, [
                np.log(2) * beta_ref[0], 0, 0], disp=False)
            beta_npp0[i, b] = 1 + np.abs(beta_npp0[i, b])
            CO2_comp[i, b] = np.abs(CO2_comp[i, b])


            def err(var):
                K2 = ((2 * CO2_cmip5[150] - CO2_comp[i, b]) - beta_npp0[i, b] * (
                            CO2_cmip5[150] - CO2_comp[i, b])) / ((beta_npp0[i, b] - 1) * (
                            CO2_cmip5[150] - CO2_comp[i, b]) * (2 * CO2_cmip5[150] -
                                                                CO2_comp[i, b]))
                K1 = (1 + K2 * (CO2_cmip5[150] - CO2_comp[i, b])) / (
                            CO2_cmip5[150] - CO2_comp[i, b])
                carb = (K1 * (CO2_all[:] - CO2_comp[i, b])) / (
                            1 + K2 * (CO2_all[:] - CO2_comp[i, b]))
                clim = 1 + gamma_nppT[i, b] * (tas_all[:, i] - np.mean(tas_ctrl[:, i])) + \
                       var[0] * (pr_all[:, i] - np.mean(pr_ctrl[:, i]))
                return np.sum((ratio - carb * clim) ** 2)


            [gamma_nppP[i, b]] = fmin(err, [0], disp=False)

        # for RH rate
        ratio = (RH_all / CSOIL0_all)[:, i, b] / np.mean((RH_ctrl / CSOIL0_ctrl)[:, i, b])
        ratio_dec = (RH_dec / CSOIL0_dec)[:, i, b] / np.mean(
            (RH_ctrl / CSOIL0_ctrl)[:, i, b])
        if (mod_LSNKrho == 'exp'):
            def err(var):
                carb = 1 + np.abs(var[0]) * (
                            FINPUT_dec[:, i, b] - np.mean(FINPUT_ctrl[:, i, b]))
                clim = np.exp(var[1] * (tas_dec[:, i] - np.mean(tas_ctrl[:, i])))
                return np.sum((ratio_dec - carb * clim) ** 2)


            first_guess = (ratio_dec[-131] - 1) / (
                        FINPUT_dec[-131, i, b] - np.mean(FINPUT_ctrl[:, i, b]))
            [prim_rho[i, b], gamma_rhoT[i, b]] = fmin(err, [first_guess,
                                                            np.log(Q10_ref[0]) / 10],
                                                      disp=False)
            prim_rho[i, b] = np.abs(prim_rho[i, b])


            def err(var):
                carb = 1 + prim_rho[i, b] * (
                            FINPUT_all[:, i, b] - np.mean(FINPUT_ctrl[:, i, b]))
                clim = np.exp(
                    gamma_rhoT[i, b] * (tas_all[:, i] - np.mean(tas_ctrl[:, i]))) * (
                                   1 + var[0] * (pr_all[:, i] - np.mean(pr_ctrl[:, i])))
                return np.sum((ratio - carb * clim) ** 2)


            [gamma_rhoP[i, b]] = fmin(err, [0], disp=False)
        elif (mod_LSNKrho == 'gauss'):
            def err(var):
                carb = 1 + np.abs(var[0]) * (
                            FINPUT_dec[:, i, b] - np.mean(FINPUT_ctrl[:, i, b]))
                clim = np.exp(
                    np.abs(var[1]) * (tas_dec[:, i] - np.mean(tas_ctrl[:, i])) - np.abs(
                        var[2]) * (tas_dec[:, i] - np.mean(tas_ctrl[:, i])) ** 2)
                return np.sum((ratio_dec - carb * clim) ** 2)


            first_guess = (ratio_dec[-131] - 1) / (
                        FINPUT_dec[-131, i, b] - np.mean(FINPUT_ctrl[:, i, b]))
            [prim_rho[i, b], gamma_rhoT1[i, b], gamma_rhoT2[i, b]] = fmin(err,
                                                                          [first_guess,
                                                                           np.log(Q10_ref[
                                                                                      0]) / 10,
                                                                           0], disp=False)
            prim_rho[i, b] = np.abs(prim_rho[i, b])
            gamma_rhoT1[i, b] = np.abs(gamma_rhoT1[i, b])
            gamma_rhoT2[i, b] = -np.abs(gamma_rhoT2[i, b])


            def err(var):
                carb = 1 + prim_rho[i, b] * (
                            FINPUT_all[:, i, b] - np.mean(FINPUT_ctrl[:, i, b]))
                clim = np.exp(
                    gamma_rhoT1[i, b] * (tas_all[:, i] - np.mean(tas_ctrl[:, i])) +
                    gamma_rhoT2[i, b] * (
                                tas_all[:, i] - np.mean(tas_ctrl[:, i])) ** 2) * (
                                   1 + var[0] * (pr_all[:, i] - np.mean(pr_ctrl[:, i])))
                return np.sum((ratio - carb * clim) ** 2)


            [gamma_rhoP[i, b]] = fmin(err, [0], disp=False)

# arbitrary values of parameters
# avoid diverging value of sigma
if (mod_LSNKnpp == 'hyp'):
    beta_npp0[beta_npp0 == 1] = 1.000001
# urban: no change to preindustrial
if (nb_biome > 1) & ((mod_biomeURB == 'URB') | mod_biomeV3):
    if (mod_LSNKnpp == 'log'):
        beta_npp[:, biome_index["urb"]] = 0
    elif (mod_LSNKnpp == 'hyp'):
        beta_npp0[:, biome_index["urb"]] = 1.000001
        CO2_comp[:, biome_index["urb"]] = 1E18
    gamma_nppT[:, biome_index["urb"]] = 0
    gamma_nppP[:, biome_index["urb"]] = 0
    prim_rho[:, biome_index["urb"]] = 0
    if (mod_LSNKrho == 'exp'):
        gamma_rhoT[:, biome_index["urb"]] = 0
        gamma_rhoP[:, biome_index["urb"]] = 0
    elif (mod_LSNKrho == 'gauss'):
        gamma_rhoT1[:, biome_index["urb"]] = 0
        gamma_rhoT2[:, biome_index["urb"]] = 0
        gamma_rhoP[:, biome_index["urb"]] = 0

# load pre-processed CMIP5 results for specified model
# data related to wildfires
for sim in ['ctrl', 'upct', 'fxcl', 'fdbk']:
    for VAR in ['EFIRE', 'CBURNT']:
        exec(VAR + '_' + sim + ' = np.zeros([140,nb_regionI,len(bio)], dtype=dty)')
        if (mod_EFIREtrans != ''):
            for b in range(len(bio)):
                if os.path.isfile(
                        'data/BioBurn_CMIP5/#DATA.BioBurn_' + mod_EFIREtrans + '.140yr_114reg1.' + sim + '_' + VAR + '_' +
                        bio[b] + '.csv'):
                    TMP = np.array([line for line in csv.reader(open(
                        'data/BioBurn_CMIP5/#DATA.BioBurn_' + mod_EFIREtrans + '.140yr_114reg1.' + sim + '_' + VAR + '_' +
                        bio[b] + '.csv', 'r'))], dtype=dty)
                    for i in range(1, 114 + 1):
                        exec(VAR + '_' + sim + '[:,regionI_index[i],b] += TMP[:,i-1]')
    for var in ['tas2', 'pr2']:
        exec(var + '_' + sim + ' = np.zeros([140,nb_regionI], dtype=dty)')
        if (mod_EFIREtrans != ''):
            TMP = np.array([line for line in csv.reader(open(
                'data/BioBurn_CMIP5/#DATA.BioBurn_' + mod_EFIREtrans + '.140yr_114reg1.' + sim + '_' + var + '.csv',
                'r'))], dtype=dty)
            TMP2 = np.array([line for line in csv.reader(open(
                'data/BioBurn_CMIP5/#DATA.BioBurn_' + mod_EFIREtrans + '.140yr_114reg1.AREA2.csv',
                'r'))], dtype=dty)
            for i in range(1, 114 + 1):
                exec(var + '_' + sim + '[:,regionI_index[i]] += TMP[:,i-1]*TMP2[:,i-1]')
            TMP = np.zeros([140, nb_regionI], dtype=dty)
            for i in range(1, 114 + 1):
                exec('TMP[:,regionI_index[i]] += TMP2[:,i-1]')
            exec(var + '_' + sim + ' /= TMP')
            exec(
                var + '_' + sim + '[np.isnan(' + var + '_' + sim + ')|np.isinf(' + var + '_' + sim + ')] = 0')

# complete with arbitrary values
# if no shrubs in the model
if (nb_biome > 1) & (np.sum(CBURNT_ctrl[:, :, bio.index("shr")]) == 0) & (
        mod_biomeSHR == 'SHR'):
    for sim in ['ctrl', 'upct', 'fxcl', 'fdbk']:
        exec(
            'CBURNT_' + sim + '[:,:,bio.index("shr")] = CBURNT_' + sim + '[:,:,bio.index("gra")]')
        for VAR in ['EFIRE']:
            exec(
                VAR + '_' + sim + '[:,:,bio.index("shr")] = 0.85*' + VAR + '_' + sim + '[:,:,bio.index("gra")] + 0.15*' + VAR + '_' + sim + '[:,:,bio.index("for")]*CBURNT_' + sim + '[:,:,bio.index("gra")]/CBURNT_' + sim + '[:,:,bio.index("for")]')
            exec('TMP = (CBURNT_' + sim + '[:,:,bio.index("for")] == 0)')
            exec(
                VAR + '_' + sim + '[:,:,bio.index("shr")][TMP] = ' + VAR + '_' + sim + '[:,:,bio.index("gra")][TMP]')
# if no crops in the model
if (nb_biome > 1) & (np.sum(CBURNT_ctrl[:, :, bio.index("cro")]) == 0):
    for sim in ['ctrl', 'upct', 'fxcl', 'fdbk']:
        for VAR in ['CBURNT', 'EFIRE']:
            exec(
                VAR + '_' + sim + '[:,:,bio.index("cro")] = ' + VAR + '_' + sim + '[:,:,bio.index("gra")]')
# if no pastures in the model
if (nb_biome > 1) & (np.sum(CBURNT_ctrl[:, :, bio.index("pas")]) == 0):
    for sim in ['ctrl', 'upct', 'fxcl', 'fdbk']:
        exec(
            'CBURNT_' + sim + '[:,:,bio.index("pas")] = CBURNT_' + sim + '[:,:,bio.index("gra")]')
        for VAR in ['EFIRE']:
            exec(
                VAR + '_' + sim + '[:,:,bio.index("pas")] = 0.60*' + VAR + '_' + sim + '[:,:,bio.index("gra")] + 0.40*' + VAR + '_' + sim + '[:,:,bio.index("des")]*CBURNT_' + sim + '[:,:,bio.index("gra")]/CBURNT_' + sim + '[:,:,bio.index("des")]')
            exec('TMP = (CBURNT_' + sim + '[:,:,bio.index("des")] == 0)')
            exec(
                VAR + '_' + sim + '[:,:,bio.index("pas")][TMP] = ' + VAR + '_' + sim + '[:,:,bio.index("gra")][TMP]')
# if no urban in the model
if (nb_biome > 1) & (np.sum(CBURNT_ctrl[:, :, bio.index("urb")]) == 0) & (
        (mod_biomeURB == 'URB') | mod_biomeV3):
    for sim in ['ctrl', 'upct', 'fxcl', 'fdbk']:
        for VAR in ['CBURNT', 'EFIRE']:
            exec(
                VAR + '_' + sim + '[:,:,bio.index("urb")] = ' + VAR + '_' + sim + '[:,:,bio.index("des")]')

# aggregate biomes
for sim in ['ctrl', 'upct', 'fxcl', 'fdbk']:
    for VAR in ['EFIRE', 'CBURNT']:
        exec('TMP = ' + VAR + '_' + sim + '.copy()')
        exec(VAR + '_' + sim + ' = np.zeros([140,nb_regionI,nb_biome], dtype=dty)')
        for b in range(len(bio)):
            exec(VAR + '_' + sim + '[:,:,biome_index[bio[b]]] += TMP[:,:,b]')
# aggregate all experiments
CO2_all = np.array(list(CO2_cmip5[150] * 1.01 ** np.arange(140)) + list(
    CO2_cmip5[150] * 1.01 ** np.arange(140)) + list(CO2_cmip5[150] * np.ones(140)),
                   dtype=dty)
for VAR in ['EFIRE', 'CBURNT'] + ['tas2', 'pr2']:
    exec(
        VAR + '_all = np.array(list(' + VAR + '_upct)+list(' + VAR + '_fxcl)+list(' + VAR + '_fdbk), dtype=dty)')
# decadal means
for VAR in ['EFIRE', 'CBURNT'] + ['tas2', 'pr2']:
    exec(
        VAR + '_dec= np.zeros([3*(140-10)]+list(np.shape(' + VAR + '_all)[1:]), dtype=dty)')
    for t in range(130):
        exec(VAR + '_dec[t,...] = np.mean(' + VAR + '_all[t:t+10,...],0)')
        exec(VAR + '_dec[t+140-10,...] = np.mean(' + VAR + '_all[140+t:140+t+10,...],0)')
        exec(
            VAR + '_dec[t+2*140-2*10,...] = np.mean(' + VAR + '_all[2*140+t:2*140+t+10,...],0)')

# definition of parameters
# sensitivity of ignition rate to climate {/ppm}&{/K}&{/mm}
gamma_igniC = np.zeros([nb_regionI, nb_biome], dtype=dty)
gamma_igniT = np.zeros([nb_regionI, nb_biome], dtype=dty)
gamma_igniP = np.zeros([nb_regionI, nb_biome], dtype=dty)

# fit of parameters
for i in range(1, nb_regionI):
    for b in range(nb_biome):

        # for ignition rate
        if (mod_EFIREtrans != ''):
            ratio = (EFIRE_all / CBURNT_all)[:, i, b] / np.mean(
                (EFIRE_ctrl / CBURNT_ctrl)[:, i, b])
            ratio_dec = (EFIRE_dec / CBURNT_dec)[:, i, b] / np.mean(
                (EFIRE_ctrl / CBURNT_ctrl)[:, i, b])


            def err(var):
                carb = 1
                clim = 1 + var[0] * (CO2_dec[:] - CO2_cmip5[150]) + np.abs(var[1]) * (
                            tas2_dec[:, i] - np.mean(tas2_ctrl[:, i])) - np.abs(
                    var[2]) * (pr2_dec[:, i] - np.mean(pr2_ctrl[:, i]))
                return np.sum((ratio_dec - carb * clim) ** 2)


            first_guess = (ratio_dec[-131] - 1) / (CO2_dec[-131] - CO2_cmip5[150])
            [gamma_igniC[i, b], gamma_igniT[i, b], gamma_igniP[i, b]] = fmin(err,
                                                                             [first_guess,
                                                                              0, 0],
                                                                             disp=False)
            gamma_igniT[i, b] = np.abs(gamma_igniT[i, b])
            gamma_igniP[i, b] = -np.abs(gamma_igniP[i, b])

# arbitrary values of parameters
# crops: no wildfire
if (nb_biome > 1):
    gamma_igniC[:, biome_index["cro"]] = 0.
    gamma_igniT[:, biome_index["cro"]] = 0.
    gamma_igniP[:, biome_index["cro"]] = 0.
# urban: no wildfire
if (nb_biome > 1) & ((mod_biomeURB == 'URB') | mod_biomeV3):
    gamma_igniC[:, biome_index["urb"]] = 0.
    gamma_igniT[:, biome_index["urb"]] = 0.
    gamma_igniP[:, biome_index["urb"]] = 0.

# -----------------
# 1.3.5. Permafrost
# -----------------

# permafrost regions
regionPF = ['EUAS', 'NAM']
regionPF_name = ['Eurasia', 'North America']
nb_regionPF = len(regionPF)


# function for speed of thawing/freezing
def f_v_PF(pthaw_bar, pthaw):
    v = v_thaw * (pthaw_bar - pthaw > 0) + v_froz * (pthaw_bar - pthaw < 0)
    return np.array(v, dtype=dty)


# load other calibrated parameters
if (mod_EPFmain == 'JSBACH'):

    # regional weighting for local temperature {.}
    w_reg_lstPF = np.array([1.8568359, 1.9470703], dtype=dty)

    # sensitivities of local heterotrophic respiration {/K}&{/K2}&{.}
    gamma_rhoPF1 = np.array([0.09690543, 0.10076028], dtype=dty)
    gamma_rhoPF2 = -np.array([0.00206356, 0.00214174], dtype=dty)
    w_rhoPF = np.array([0.93812785, 0.73454044], dtype=dty)

    # parameters for thawing function {.}&{/K}&{.}
    pthaw_min = np.array([0.25319850, 0.19503177], dtype=dty)
    k_pthaw = np.array([2.3012887, 1.4174312], dtype=dty)
    gamma_pthaw = np.array([0.17582392, 0.28638024], dtype=dty)

    # speed of thawing/freezing {/yr}
    v_thaw = np.array([4.8940873, 0.49978620], dtype=dty)
    v_froz = np.array([0., 0.], dtype=dty)

    # weights and turnover times for the three emitting pools {.}&{yr}
    p_PF1 = np.array([0.06976970, 0.09987357], dtype=dty)
    p_PF2 = np.array([0.26396975, 0.23425585], dtype=dty)
    p_PF3 = np.array([0.66626055, 0.66587058], dtype=dty)
    tau_PF1 = np.array([7.9616249, 11.632080], dtype=dty)
    tau_PF2 = np.array([88.489543, 103.06668], dtype=dty)
    tau_PF3 = np.array([1355.6291, 1240.6807], dtype=dty)

    # initial pool of frozen carbon {GtC}
    CFROZ_0 = np.array([413.553, 277.483], dtype=dty)

elif (mod_EPFmain == 'ORCHIDEE-MICT'):

    # regional weighting for local temperature {.}
    w_reg_lstPF = np.array([1.8678710, 1.7828125], dtype=dty)

    # sensitivities of local heterotrophic respiration {/K}&{/K2}&{.}
    gamma_rhoPF1 = np.array([0.11400756, 0.10956953], dtype=dty)
    gamma_rhoPF2 = -np.array([0.00366366, 0.00345403], dtype=dty)
    w_rhoPF = np.array([3.3805183, 3.7076802], dtype=dty)

    # parameters for thawing function {.}&{/K}&{.}
    pthaw_min = np.array([0.08065349, 0.87546773], dtype=dty)
    k_pthaw = np.array([0.35689219, 8.4422619], dtype=dty)
    gamma_pthaw = np.array([1.6183032, 0.09725753], dtype=dty)

    # speed of thawing/freezing {/yr}
    v_thaw = np.array([0.28696367, 0.49953314], dtype=dty)
    v_froz = np.array([0.05289277, 0.04902789], dtype=dty)

    # weights and turnover times for the three emitting pools {.}&{yr}
    p_PF1 = np.array([0.24389266, 0.28834940], dtype=dty)
    p_PF2 = np.array([0.75610734, 0.71165060], dtype=dty)
    p_PF3 = np.array([0., 0.], dtype=dty)
    tau_PF1 = np.array([1332.6992, 2950.1099], dtype=dty)
    tau_PF2 = np.array([14953.331, 23437.897], dtype=dty)
    tau_PF3 = np.array([np.inf, np.inf], dtype=dty)

    # initial pool of frozen carbon {GtC}
    CFROZ_0 = np.array([270.56425, 118.19786], dtype=dty)

elif (mod_EPFmain == 'JULES-DeepResp'):

    # regional weighting for local temperature {.}
    w_reg_lstPF = np.array([1.9603515, 2.0255859], dtype=dty)

    # sensitivities of local heterotrophic respiration {/K}&{/K2}&{.}
    gamma_rhoPF1 = np.array([0.16817739, 0.14376888], dtype=dty)
    gamma_rhoPF2 = -np.array([0.00531266, 0.00379630], dtype=dty)
    w_rhoPF = np.array([1.3748762, 1.5823555], dtype=dty)

    # parameters for thawing function {.}&{/K}&{.}
    pthaw_min = np.array([0.77246763, 0.95893920], dtype=dty)
    k_pthaw = np.array([2.0529016, 1.9143394], dtype=dty)
    gamma_pthaw = np.array([0.1530148, 0.14457936], dtype=dty)

    # speed of thawing/freezing {/yr}
    v_thaw = np.array([0.26619401, 0.59899092], dtype=dty)
    v_froz = np.array([0.04102776, 0.03296677], dtype=dty)

    # weights and turnover times for the three emitting pools {.}&{yr}
    p_PF1 = np.array([0.03366769, 0.65114842], dtype=dty)
    p_PF2 = np.array([0.96633231, 0.34885158], dtype=dty)
    p_PF3 = np.array([0., 0.], dtype=dty)
    tau_PF1 = np.array([156.38595, 1969.9072], dtype=dty)
    tau_PF2 = np.array([3157.9293, 5370.1557], dtype=dty)
    tau_PF3 = np.array([np.inf, np.inf], dtype=dty)

    # initial pool of frozen carbon {GtC}
    CFROZ_0 = np.array([480.92224, 175.53549], dtype=dty)

elif (mod_EPFmain == 'JULES-SuppressResp'):

    # regional weighting for local temperature {.}
    w_reg_lstPF = np.array([1.9590820, 2.0242187], dtype=dty)

    # sensitivities of local heterotrophic respiration {/K}&{/K2}&{.}
    gamma_rhoPF1 = np.array([0.13472092, 0.11193608], dtype=dty)
    gamma_rhoPF2 = -np.array([0.00498800, 0.00339437], dtype=dty)
    w_rhoPF = np.array([0.66180279, 2.1853523], dtype=dty)

    # parameters for thawing function {.}&{/K}&{.}
    pthaw_min = np.array([0.87020023, 1.1649739], dtype=dty)
    k_pthaw = np.array([1.6305246, 1.6836054], dtype=dty)
    gamma_pthaw = np.array([0.17978647, 0.19005963], dtype=dty)

    # speed of thawing/freezing {/yr}
    v_thaw = np.array([0.37339128, 0.49223681], dtype=dty)
    v_froz = np.array([0.06128438, 0.04292502], dtype=dty)

    # weights and turnover times for the three emitting pools {.}&{yr}
    p_PF1 = np.array([1., 1.], dtype=dty)
    p_PF2 = np.array([0., 0.], dtype=dty)
    p_PF3 = np.array([0., 0.], dtype=dty)
    tau_PF1 = np.array([9433.4705, 4322.9007], dtype=dty)
    tau_PF2 = np.array([np.inf, np.inf], dtype=dty)
    tau_PF3 = np.array([np.inf, np.inf], dtype=dty)

    # initial pool of frozen carbon {GtC}
    CFROZ_0 = np.array([372.75814, 121.20251], dtype=dty)

elif (mod_EPFmain == ''):

    # regional weighting for local temperature {.}
    w_reg_lstPF = np.array([1., 1.], dtype=dty)

    # sensitivities of local heterotrophic respiration {/K}&{/K2}&{.}
    gamma_rhoPF1 = np.array([0., 0.], dtype=dty)
    gamma_rhoPF2 = -np.array([0., 0.], dtype=dty)
    w_rhoPF = np.array([1., 1.], dtype=dty)

    # parameters for thawing function {.}&{/K}&{.}
    pthaw_min = np.array([0., 0.], dtype=dty)
    k_pthaw = np.array([1., 1.], dtype=dty)
    gamma_pthaw = np.array([0., 0.], dtype=dty)

    # speed of thawing/freezing {/yr}
    v_thaw = np.array([1., 1.], dtype=dty)
    v_froz = np.array([0., 0.], dtype=dty)

    # weights and turnover times for the three emitting pools {.}&{yr}
    p_PF1 = np.array([0., 0.], dtype=dty)
    p_PF2 = np.array([0., 0.], dtype=dty)
    p_PF3 = np.array([0., 0.], dtype=dty)
    tau_PF1 = np.array([np.inf, np.inf], dtype=dty)
    tau_PF2 = np.array([np.inf, np.inf], dtype=dty)
    tau_PF3 = np.array([np.inf, np.inf], dtype=dty)

    # initial pool of frozen carbon {GtC}
    CFROZ_0 = np.array([0., 0.], dtype=dty)

# =============
# 1.4. LAND-USE
# =============

# -------------
# 1.4.1. TRENDY
# -------------

# basic biomes of aggregation
bio = ['des', 'for', 'shr', 'gra', 'cro', 'pas', 'urb']

# load pre-processed TRENDY results for specified model
# data related to preindustrial C-cycle
for VAR in ['AGB', 'BGB']:
    TMP = np.array([line for line in csv.reader(open(
        'data/WoodUse_TRENDYv2/#DATA.WoodUse_' + mod_ELUCagb + '.1910s_114reg1_7bio.' + VAR + '.csv',
        'r'))], dtype=dty)
    exec(VAR + '_pi = np.zeros([nb_regionI,len(bio)], dtype=dty)')
    for i in range(1, 114 + 1):
        exec(VAR + '_pi[regionI_index[i],:] += TMP[i-1,:]')

# complete with arbitrary values
# if no shrubs in the model
if (nb_biome > 1) & (np.sum((AGB_pi + BGB_pi)[:, bio.index("shr")]) == 0) & (
        mod_biomeSHR == 'SHR'):
    for VAR in ['AGB', 'BGB']:
        exec(
            VAR + '_pi[:,bio.index("shr")] = 0.85*' + VAR + '_pi[:,bio.index("gra")] + 0.15*' + VAR + '_pi[:,bio.index("for")]')
# if no crops in the model
if (nb_biome > 1) & (np.sum((AGB_pi + BGB_pi)[:, bio.index("cro")]) == 0):
    for VAR in ['AGB', 'BGB']:
        exec(VAR + '_pi[:,bio.index("cro")] = ' + VAR + '_pi[:,bio.index("gra")]')
# if no pastures in the model
if (nb_biome > 1) & (np.sum((AGB_pi + BGB_pi)[:, bio.index("pas")]) == 0):
    for VAR in ['AGB', 'BGB']:
        exec(
            VAR + '_pi[:,bio.index("pas")] = 0.60*' + VAR + '_pi[:,bio.index("gra")] + 0.40*' + VAR + '_pi[:,bio.index("des")]')
# if no urban in the model
if (nb_biome > 1) & (np.sum((AGB_pi + BGB_pi)[:, bio.index("urb")]) == 0) & (
        (mod_biomeURB == 'URB') | mod_biomeV3):
    for VAR in ['AGB', 'BGB']:
        exec(VAR + '_pi[:,bio.index("urb")] = ' + VAR + '_pi[:,bio.index("des")]')

# aggregate biomes
for VAR in ['AGB', 'BGB']:
    exec('TMP = ' + VAR + '_pi.copy()')
    exec(VAR + '_pi = np.zeros([nb_regionI,nb_biome], dtype=dty)')
    for b in range(len(bio)):
        exec(VAR + '_pi[:,biome_index[bio[b]]] += TMP[:,b]')

# calculate preindustrial flux parameters
p_AGB = AGB_pi / (AGB_pi + BGB_pi)
p_AGB[np.isnan(p_AGB)] = 0

# ---------------
# 1.4.2. Wood-Use
# ---------------

# characteristic time for shifting cultivation {yr}
# adjusted; 15yr is from [Hurtt et al., 2006]
tau_shift = np.array([15], dtype=dty)

# aggregation of fractions for wood-use {.}
# (n=0) for slash; (n=1) for burnt; (n=2,3) for wood decay
# base data from [Earles et al., 2012]
for n in range(4):
    exec('p_HWP' + str(n) + ' = np.zeros([nb_regionI,nb_biome], dtype=dty)')
TMP = np.array([line for line in csv.reader(
    open('data/WoodUse_Earles/#DATA.WoodUse_Earles.2000s_114reg1_(4use).HWP.csv', 'r'))][
               1:], dtype=dty)
for i in range(1, 114 + 1):
    for n in range(4):
        exec('p_HWP' + str(n) + '[regionI_index[i],:] += TMP[i-1,n]')
# biomass-burning of non-merchantable AGB
if (mod_EHWPbb == 'high'):
    p_HWP1 += 0.5 * p_HWP0
    p_HWP0 -= 0.5 * p_HWP0
elif (mod_EHWPbb == 'low'):
    p_HWP1 += 0.0 * p_HWP0
    p_HWP0 -= 0.0 * p_HWP0
# remove some products for some biomes
if (nb_biome > 1):
    for bio in ['des', 'shr', 'gra']:
        p_HWP2[:, biome_index[bio]] = p_HWP3[:, biome_index[bio]] = 0
    for bio in ['cro', 'pas']:
        p_HWP1[:, biome_index[bio]] = p_HWP2[:, biome_index[bio]] = p_HWP3[:,
                                                                    biome_index[bio]] = 0
    if (mod_biomeURB == 'URB') | mod_biomeV3:
        p_HWP1[:, biome_index["urb"]] = p_HWP2[:, biome_index["urb"]] = p_HWP3[:,
                                                                        biome_index[
                                                                            "urb"]] = 0
# normalize
TMP = p_HWP0 + p_HWP1 + p_HWP2 + p_HWP3
for n in range(4):
    exec('p_HWP' + str(n) + ' /= TMP')
    exec('p_HWP' + str(n) + '[np.isinf(p_HWP' + str(n) + ')|np.isnan(p_HWP' + str(
        n) + ')] = 0')
p_HWP0 = 1 - p_HWP1 - p_HWP2 - p_HWP3

# fraction of actually burnt HWP1 {}
# used to adjust non-CO2 anthropogenic BB emissions
p_HWP1_bb = 0.5

# wood-use decay constants {yr}
# based roughly on [Houghton et Hackler, 2001]
if (mod_EHWPtau == 'Houghton2001'):
    tau_HWP1 = np.array([1], dtype=dty)
    tau_HWP2 = np.array([10], dtype=dty)
    tau_HWP3 = np.array([100], dtype=dty)
# based on [Earles et al., 2012]
elif (mod_EHWPtau == 'Earles2012'):
    tau_HWP1 = np.array([0.5], dtype=dty)
    tau_HWP2 = np.array([2], dtype=dty)
    tau_HWP3 = np.array([30], dtype=dty)

## decay times variations
## fast: 20% remaining after 0.8*tau
if mod_EHWPspeed == 'fast':
    tau_HWP1 *= 0.8 / -np.log(0.2)
    tau_HWP2 *= 0.8 / -np.log(0.2)
    tau_HWP3 *= 0.8 / -np.log(0.2)
## fast: 30% remaining after 1.5*tau
elif mod_EHWPspeed == 'slow':
    tau_HWP1 *= 1.5 / -np.log(0.3)
    tau_HWP2 *= 1.5 / -np.log(0.3)
    tau_HWP3 *= 1.5 / -np.log(0.3)

# IRF for wood product pools {.}
## forced to be exponential in this version
mod_EHWPfct = 'exp'
if (mod_EHWPfct == 'gamma'):
    def f_HWP(tau):
        err = lambda a0: gammainc(a0, (a0 + 1) / 2.5) - 0.05
        a0 = fsolve(err, 0.1)
        HWP = gammainc(a0, tau * (a0 + 1) / np.arange(ind_final + 1))
        HWP[0] = 1.
        return np.array(HWP, dtype=dty)
elif (mod_EHWPfct == 'lin'):
    def f_HWP(tau):
        HWP = (1 - np.arange(ind_final + 1) / tau) * (np.arange(ind_final + 1) < tau)
        return np.array(HWP, dtype=dty)
elif (mod_EHWPfct == 'exp'):
    def f_HWP(tau):
        HWP = np.exp(-np.arange(ind_final + 1) / tau)
        return np.array(HWP, dtype=dty)

# IRF for wood product fluxes ratio {/yr}
for N in ['1', '2', '3']:
    exec('r_HWP' + N + ' = np.zeros([ind_final+1], dtype=dty)')
    exec(
        'r_HWP' + N + '[1:] = -(f_HWP(tau_HWP' + N + ')[1:]-f_HWP(tau_HWP' + N + ')[:-1])/(f_HWP(tau_HWP' + N + ')[1:]+f_HWP(tau_HWP' + N + ')[:-1])/0.5')
    exec('r_HWP' + N + '[np.isnan(r_HWP' + N + ')|np.isinf(r_HWP' + N + ')] = 1')
    # exec('r_HWP'+N+'[r_HWP'+N+'>1] = 1')
    exec('r_HWP' + N + '[0] = 0')

# -----------
# 1.4.3. GFED
# -----------

# biomes of GFED and correspondence
bio = ['def', 'for', 'woo', 'sav', 'agr', 'pea']
BIO = {'def': 'for', 'for': 'for', 'woo': 'shr', 'sav': 'gra', 'agr': 'cro', 'pea': 'pea'}

# load GFED v3.1 BB emissions {TgX/GtC}&{Tg/GtC}
# from [Randerson et al., 2013]
for VAR in ['CO2'] + ['CH4', 'NOX', 'CO', 'VOC', 'N2O', 'SO2', 'NH3', 'BC', 'OC']:
    exec('alpha_BB_' + VAR + ' = np.zeros([nb_regionI,nb_biome], dtype=dty)')
    for b in range(len(bio) - 1):
        if os.path.isfile(
                'data/BioBurn_GFED/#DATA.BioBurn_GFED.1997-2011_114reg1.E' + VAR + '_' +
                bio[b] + '.csv'):
            TMP = np.array([line for line in csv.reader(open(
                'data/BioBurn_GFED/#DATA.BioBurn_GFED.1997-2011_114reg1.E' + VAR + '_' +
                bio[b] + '.csv', 'r'))], dtype=dty)
            for i in range(1, 114 + 1):
                exec(
                    'alpha_BB_' + VAR + '[regionI_index[i],biome_index[BIO[bio[b]]]] += np.sum(TMP[:,i-1])')
                # set pastures and baresoil as grasslands
                if (bio[b] == 'sav'):
                    exec(
                        'alpha_BB_' + VAR + '[regionI_index[i],biome_index["des"]] += np.sum(TMP[:,i-1])')
                    exec(
                        'alpha_BB_' + VAR + '[regionI_index[i],biome_index["pas"]] += np.sum(TMP[:,i-1])')
for VAR in ['CH4', 'NOX', 'CO', 'VOC', 'N2O', 'SO2', 'NH3', 'BC', 'OC'] + ['CO2']:
    exec('alpha_BB_' + VAR + ' /= alpha_BB_CO2')
    exec(
        'alpha_BB_' + VAR + '[np.isnan(alpha_BB_' + VAR + ')|np.isinf(alpha_BB_' + VAR + ')] = 0')
    for b in range(nb_biome):
        exec(
            'alpha_BB_' + VAR + '[:,b][alpha_BB_' + VAR + '[:,b]==0] = np.sum(alpha_BB_' + VAR + '[:,b])/np.sum(alpha_BB_' + VAR + '[:,b]!=0)')
    exec(
        'alpha_BB_' + VAR + '[np.isnan(alpha_BB_' + VAR + ')|np.isinf(alpha_BB_' + VAR + ')] = 0')

##################################################
#   2. METHANE
##################################################

# ===============
# 2.1. ATMOSPHERE
# ===============

# conversion of CH4 from {ppb} to {TgC}
alpha_CH4 = 0.1765 * np.array([12.0], dtype=dty)

# historic CH4 from IPCC-AR5 {ppb}
# from [IPCC WG1, 2013] annexe 2
CH4_ipcc = np.ones([311 + 1], dtype=dty) * np.nan
TMP = np.array([line for line in csv.reader(
    open('data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1750-2011.CH4.csv', 'r'))],
               dtype=dty)
CH4_ipcc[50:] = TMP[:, 0]
CH4_0 = np.array([CH4_ipcc[50]], dtype=dty)

# historic CH4 from CMIP5 {ppb}
# from [Meinshausen et al., 2011]
CH4_cmip5 = np.ones([305 + 1], dtype=dty) * np.nan
TMP = np.array([line for line in csv.reader(
    open('data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005.CH4.csv', 'r'))], dtype=dty)
CH4_cmip5[65:] = TMP[:, 0]

# historic CH4 from AGAGE {ppb}
# from [Prinn et al., 2013] updated from the website
CH4_agage = np.ones([313 + 1], dtype=dty) * np.nan
TMP = np.array([line for line in csv.reader(
    open('data/HistAtmo_AGAGE/#DATA.HistAtmo_AGAGE.1987-2013.CH4_global.csv', 'r'))],
               dtype=dty)
CH4_agage[287:] = TMP[:, 0]

# historic CH4 from Law Dome ice cores {ppm}
# from [Etheridge et al., 1998] and [MacFarling Meure et al., 2006]
CH4_lawdome = np.array([line for line in csv.reader(
    open('data/HistAtmo_NOAA-NCDC/#DATA.HistAtmo_NOAA-NCDC.(IceCores).CH4_lawdome.csv',
         'r'))][1:], dtype=dty)

# load RCP concentrations {ppb}
# from [Meinshausen et al., 2011]
CH4_rcp = np.ones([800 + 1, 6], dtype=dty) * np.nan
n = -1
for rcp in ['rcp26', 'rcp45', 'rcp60', 'rcp85', 'rcp45to26', 'rcp60to45']:
    n += 1
    TMP = np.array([line for line in csv.reader(
        open('data/Scenario_ECP/#DATA.Scenario_ECP.2000-2500.' + rcp + '_CH4.csv', 'r'))],
                   dtype=dty)
    CH4_rcp[300:, n] = TMP[:, 0]

# ==============
# 2.2. CHEMISTRY
# ==============

# ----------------
# 2.2.1. Lifetimes
# ----------------

# preindustrial lifetime of CH4 {yr}
# from [Prather et al., 2012]
# rescale of strato lifetime from [Prather et al., 2015]
tau_CH4_hv = 1.06 * np.array([120.], dtype=dty)
tau_CH4_soil = np.array([150.], dtype=dty)
tau_CH4_ocean = np.array([200.], dtype=dty)
# average from [Prather et al., 2012]
if (mod_OHSNKtau == 'Prather2012'):
    tau_CH4_OH = np.array([11.2], dtype=dty)
# variations from [Naik et al., 2013] (table 1)
elif (mod_OHSNKtau == 'CESM-CAM-superfast'):
    tau_CH4_OH = np.array([8.4], dtype=dty) * 11.2 / 9.7
elif (mod_OHSNKtau == 'CICERO-OsloCTM2'):
    tau_CH4_OH = np.array([10.0], dtype=dty) * 11.2 / 9.7
elif (mod_OHSNKtau == 'CMAM'):
    tau_CH4_OH = np.array([9.4], dtype=dty) * 11.2 / 9.7
elif (mod_OHSNKtau == 'EMAC'):
    tau_CH4_OH = np.array([9.1], dtype=dty) * 11.2 / 9.7
elif (mod_OHSNKtau == 'GEOSCCM'):
    tau_CH4_OH = np.array([9.6], dtype=dty) * 11.2 / 9.7
elif (mod_OHSNKtau == 'GDFL-AM3'):
    tau_CH4_OH = np.array([9.4], dtype=dty) * 11.2 / 9.7
elif (mod_OHSNKtau == 'GISS-E2-R'):
    tau_CH4_OH = np.array([10.6], dtype=dty) * 11.2 / 9.7
elif (mod_OHSNKtau == 'GISS-E2-R-TOMAS'):
    tau_CH4_OH = np.array([9.2], dtype=dty) * 11.2 / 9.7
elif (mod_OHSNKtau == 'HadGEM2'):
    tau_CH4_OH = np.array([11.6], dtype=dty) * 11.2 / 9.7
elif (mod_OHSNKtau == 'LMDzORINCA'):
    tau_CH4_OH = np.array([10.5], dtype=dty) * 11.2 / 9.7
elif (mod_OHSNKtau == 'MIROC-CHEM'):
    tau_CH4_OH = np.array([8.7], dtype=dty) * 11.2 / 9.7
elif (mod_OHSNKtau == 'MOCAGE'):
    tau_CH4_OH = np.array([7.1], dtype=dty) * 11.2 / 9.7
elif (mod_OHSNKtau == 'NCAR-CAM-35'):
    tau_CH4_OH = np.array([9.2], dtype=dty) * 11.2 / 9.7
elif (mod_OHSNKtau == 'STOC-HadAM3'):
    tau_CH4_OH = np.array([9.1], dtype=dty) * 11.2 / 9.7
elif (mod_OHSNKtau == 'TM5'):
    tau_CH4_OH = np.array([9.9], dtype=dty) * 11.2 / 9.7
elif (mod_OHSNKtau == 'UM-CAM'):
    tau_CH4_OH = np.array([14.0], dtype=dty) * 11.2 / 9.7

# arbitrary rescaling (down) of OH lifetime
scale_OH = 0.80
tau_CH4_OH *= scale_OH

# -----------------
# 2.2.2. Holmes2013
# -----------------

# sensitivity of OH sink to climate, ozone and methane {.}
# from [Holmes et al., 2013]
if (mod_OHSNKtrans == 'Holmes2013'):
    chi_OH_Tatm = np.array([3.0], dtype=dty)
    chi_OH_Qatm = np.array([0.32], dtype=dty)
    chi_OH_O3 = np.array([-0.55], dtype=dty)
    chi_OH_CH4 = np.array([-0.31], dtype=dty)
elif (mod_OHSNKtrans == 'UCI-CTM'):
    chi_OH_Tatm = np.array([3.9], dtype=dty)
    chi_OH_Qatm = np.array([0.32], dtype=dty)
    chi_OH_O3 = np.array([-0.66], dtype=dty)
    chi_OH_CH4 = np.array([-0.363], dtype=dty)
elif (mod_OHSNKtrans == 'Oslo-CTM3'):
    chi_OH_Tatm = np.array([2.8], dtype=dty)
    chi_OH_Qatm = np.array([0.29], dtype=dty)
    chi_OH_O3 = np.array([-0.43], dtype=dty)
    chi_OH_CH4 = np.array([-0.307], dtype=dty)
elif (mod_OHSNKtrans == 'GEOS-Chem'):
    chi_OH_Tatm = np.array([2.2], dtype=dty)
    chi_OH_Qatm = np.array([0.34], dtype=dty)
    chi_OH_O3 = np.array([-0.61], dtype=dty)
    chi_OH_CH4 = np.array([-0.274], dtype=dty)
# but also from the TAR [Ehhalt et al., 2001]
elif (mod_OHSNKtrans == 'mean-OxComp'):
    chi_OH_Tatm = np.array([0.], dtype=dty)
    chi_OH_Qatm = np.array([0.], dtype=dty)
    chi_OH_O3 = np.array([0.], dtype=dty)
    chi_OH_CH4 = np.array([-0.32], dtype=dty)

# logarithmic sensitivity of OH sink to ozone precursors {.}
# kept logarithmic; from [Holmes et al., 2013]
if (mod_OHSNKfct == 'log'):
    if (mod_OHSNKtrans == 'Holmes2013'):
        chi_OH_NOX = np.array([-0.14], dtype=dty)
        chi_OH_CO = np.array([0.06], dtype=dty)
        chi_OH_VOC = np.array([0.04], dtype=dty)
    elif (mod_OHSNKtrans == 'UCI-CTM'):
        chi_OH_NOX = np.array([-0.15], dtype=dty)
        chi_OH_CO = np.array([0.066], dtype=dty)
        chi_OH_VOC = np.array([0.04], dtype=dty)
    elif (mod_OHSNKtrans == 'Oslo-CTM3'):
        chi_OH_NOX = np.array([-0.10], dtype=dty)
        chi_OH_CO = np.array([0.05], dtype=dty)
        chi_OH_VOC = np.array([0.04], dtype=dty)
    elif (mod_OHSNKtrans == 'GEOS-Chem'):
        chi_OH_NOX = np.array([-0.16], dtype=dty)
        chi_OH_CO = np.array([0.065], dtype=dty)
        chi_OH_VOC = np.array([0.04], dtype=dty)
    elif (mod_OHSNKtrans == 'mean-OxComp'):
        chi_OH_NOX = np.array([-0.137], dtype=dty)
        chi_OH_CO = np.array([0.11], dtype=dty)
        chi_OH_VOC = np.array([0.047], dtype=dty)

    # linear sensitivity {./{TgN/yr}}&{./{TgC/yr}}&{./{Tg/yr}}
# rescaled using the TAR [Ehhalt et al., 2001]
elif (mod_OHSNKfct == 'lin'):
    if (mod_OHSNKtrans == 'Holmes2013'):
        chi_OH_NOX = np.array([0.0042 * -0.14 / -0.137], dtype=dty)
        chi_OH_CO = np.array([-1.05E-4 * 0.06 / 0.11], dtype=dty) * 28 / 12.
        chi_OH_VOC = np.array([-3.14E-4 * 0.04 / 0.047], dtype=dty)
    elif (mod_OHSNKtrans == 'UCI-CTM'):
        chi_OH_NOX = np.array([0.0042 * -0.15 / -0.137], dtype=dty)
        chi_OH_CO = np.array([-1.05E-4 * 0.066 / 0.11], dtype=dty) * 28 / 12.
        chi_OH_VOC = np.array([-3.14E-4 * 0.04 / 0.047], dtype=dty)
    elif (mod_OHSNKtrans == 'Oslo-CTM3'):
        chi_OH_NOX = np.array([0.0042 * -0.10 / -0.137], dtype=dty)
        chi_OH_CO = np.array([-1.05E-4 * 0.05 / 0.11], dtype=dty) * 28 / 12.
        chi_OH_VOC = np.array([-3.14E-4 * 0.04 / 0.047], dtype=dty)
    elif (mod_OHSNKtrans == 'GEOS-Chem'):
        chi_OH_NOX = np.array([0.0042 * -0.16 / -0.137], dtype=dty)
        chi_OH_CO = np.array([-1.05E-4 * 0.065 / 0.11], dtype=dty) * 28 / 12.
        chi_OH_VOC = np.array([-3.14E-4 * 0.04 / 0.047], dtype=dty)
    elif (mod_OHSNKtrans == 'mean-OxComp'):
        chi_OH_NOX = np.array([0.0042], dtype=dty)
        chi_OH_CO = np.array([-1.05E-4], dtype=dty) * 28 / 12.
        chi_OH_VOC = np.array([-3.14E-4], dtype=dty)

    # constant parameters for OH sink function {.}&{.}&{K}&{DU}
# from [Holmes et al., 2013] (supp.)
k_Tatm = 0.94  # +/- 0.1
k_Qatm = 1.5  # +/- 0.1
Tatm_0 = 251.  # +/- 1
# from [Jacobson, 2005] (eq. 2.62)
k_svp = 17.67
T_svp = 243.5 - 273.15
# from [Cionni et al., 2010]
O3s_0 = 280.
# from [Skeie et al., 2011] (tab. 1)
ENOX_oh = 13.
ECO_oh = 180. * 12 / 28.
EVOC_oh = 39. + 220. + 175.

# expression of OH sink function {.}
# from [Holmes et al., 2013] (supp.)
if (mod_OHSNKfct == 'log'):

    def f_kOH(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = np.exp(chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(
            1 + D_O3s / O3s_0) + chi_OH_NOX * np.log(
            1 + ENOX / ENOX_oh) + chi_OH_CO * np.log(
            1 + ECO / ECO_oh) + chi_OH_VOC * np.log(
            1 + EVOC / EVOC_oh) + chi_OH_Tatm * np.log(
            1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(
            1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))) - 1
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dCH4(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = chi_OH_CH4 / (CH4_0 + D_CH4) * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(
                1 + D_O3s / O3s_0) + chi_OH_NOX * np.log(
                1 + ENOX / ENOX_oh) + chi_OH_CO * np.log(
                1 + ECO / ECO_oh) + chi_OH_VOC * np.log(
                1 + EVOC / EVOC_oh) + chi_OH_Tatm * np.log(
                1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1)))
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dO3s(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = chi_OH_O3 / (O3s_0 + D_O3s) * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(
                1 + D_O3s / O3s_0) + chi_OH_NOX * np.log(
                1 + ENOX / ENOX_oh) + chi_OH_CO * np.log(
                1 + ECO / ECO_oh) + chi_OH_VOC * np.log(
                1 + EVOC / EVOC_oh) + chi_OH_Tatm * np.log(
                1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1)))
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dgst(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = chi_OH_Tatm * k_Qatm * k_svp * k_Tatm / (Tatm_0 + T_svp) * np.exp(
            k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) / (1 + k_Qatm * (
                    np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1)) * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(
                1 + D_O3s / O3s_0) + chi_OH_NOX * np.log(
                1 + ENOX / ENOX_oh) + chi_OH_CO * np.log(
                1 + ECO / ECO_oh) + chi_OH_VOC * np.log(
                1 + EVOC / EVOC_oh) + chi_OH_Tatm * np.log(
                1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1)))
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dENOX(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = chi_OH_NOX / (ENOX_oh + ENOX) * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(
                1 + D_O3s / O3s_0) + chi_OH_NOX * np.log(
                1 + ENOX / ENOX_oh) + chi_OH_CO * np.log(
                1 + ECO / ECO_oh) + chi_OH_VOC * np.log(
                1 + EVOC / EVOC_oh) + chi_OH_Tatm * np.log(
                1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1)))
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dECO(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = chi_OH_CO / (ECO_oh + ECO) * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(
                1 + D_O3s / O3s_0) + chi_OH_NOX * np.log(
                1 + ENOX / ENOX_oh) + chi_OH_CO * np.log(
                1 + ECO / ECO_oh) + chi_OH_VOC * np.log(
                1 + EVOC / EVOC_oh) + chi_OH_Tatm * np.log(
                1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1)))
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dEVOC(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = chi_OH_VOC / (EVOC_oh + EVOC) * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(
                1 + D_O3s / O3s_0) + chi_OH_NOX * np.log(
                1 + ENOX / ENOX_oh) + chi_OH_CO * np.log(
                1 + ECO / ECO_oh) + chi_OH_VOC * np.log(
                1 + EVOC / EVOC_oh) + chi_OH_Tatm * np.log(
                1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1)))
        return np.array(D_kOH, dtype=dty)


    def df_kOH(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        df_1 = df_kOH_dCH4(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * D_CH4
        df_2 = df_kOH_dO3s(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * D_O3s
        df_3 = df_kOH_dgst(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * D_gst
        df_4 = df_kOH_dENOX(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * ENOX
        df_5 = df_kOH_dECO(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * ECO
        df_6 = df_kOH_dEVOC(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * EVOC
        df_tot = df_1 + df_2 + df_3 + df_4 + df_5 + df_6
        return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot),
                np.nan_to_num(df_3 / df_tot), np.nan_to_num(df_4 / df_tot),
                np.nan_to_num(df_5 / df_tot), np.nan_to_num(df_6 / df_tot)]

elif (mod_OHSNKfct == 'lin'):

    def f_kOH(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = np.exp(chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(
            1 + D_O3s / O3s_0) + chi_OH_NOX * ENOX + chi_OH_CO * ECO + chi_OH_VOC * EVOC + chi_OH_Tatm * np.log(
            1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(
            1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1))) - 1
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dCH4(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = chi_OH_CH4 / (CH4_0 + D_CH4) * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(
                1 + D_O3s / O3s_0) + chi_OH_NOX * ENOX + chi_OH_CO * ECO + chi_OH_VOC * EVOC + chi_OH_Tatm * np.log(
                1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1)))
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dO3s(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = chi_OH_O3 / (O3s_0 + D_O3s) * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(
                1 + D_O3s / O3s_0) + chi_OH_NOX * ENOX + chi_OH_CO * ECO + chi_OH_VOC * EVOC + chi_OH_Tatm * np.log(
                1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1)))
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dgst(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = chi_OH_Tatm * k_Qatm * k_svp * k_Tatm / (Tatm_0 + T_svp) * np.exp(
            k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) / (1 + k_Qatm * (
                    np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1)) * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(
                1 + D_O3s / O3s_0) + chi_OH_NOX * ENOX + chi_OH_CO * ECO + chi_OH_VOC * EVOC + chi_OH_Tatm * np.log(
                1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1)))
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dENOX(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = chi_OH_NOX * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(
                1 + D_O3s / O3s_0) + chi_OH_NOX * ENOX + chi_OH_CO * ECO + chi_OH_VOC * EVOC + chi_OH_Tatm * np.log(
                1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1)))
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dECO(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = chi_OH_CO * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(
                1 + D_O3s / O3s_0) + chi_OH_NOX * ENOX + chi_OH_CO * ECO + chi_OH_VOC * EVOC + chi_OH_Tatm * np.log(
                1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1)))
        return np.array(D_kOH, dtype=dty)


    def df_kOH_dEVOC(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        D_kOH = chi_OH_VOC * np.exp(
            chi_OH_CH4 * np.log(1 + D_CH4 / CH4_0) + chi_OH_O3 * np.log(
                1 + D_O3s / O3s_0) + chi_OH_NOX * ENOX + chi_OH_CO * ECO + chi_OH_VOC * EVOC + chi_OH_Tatm * np.log(
                1 + k_Tatm * D_gst / Tatm_0) + chi_OH_Qatm * np.log(
                1 + k_Qatm * (np.exp(k_svp * k_Tatm * D_gst / (Tatm_0 + T_svp)) - 1)))
        return np.array(D_kOH, dtype=dty)


    def df_kOH(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC):
        df_1 = df_kOH_dCH4(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * D_CH4
        df_2 = df_kOH_dO3s(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * D_O3s
        df_3 = df_kOH_dgst(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * D_gst
        df_4 = df_kOH_dENOX(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * ENOX
        df_5 = df_kOH_dECO(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * ECO
        df_6 = df_kOH_dEVOC(D_CH4, D_O3s, D_gst, ENOX, ECO, EVOC) * EVOC
        df_tot = df_1 + df_2 + df_3 + df_4 + df_5 + df_6
        return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot),
                np.nan_to_num(df_3 / df_tot), np.nan_to_num(df_4 / df_tot),
                np.nan_to_num(df_5 / df_tot), np.nan_to_num(df_6 / df_tot)]

# =========
# 2.3. LAND
# =========

# ---------------
# 2.3.1. WETCHIMP
# ---------------

# basic biomes of aggregation
bio = ['des', 'for', 'shr', 'gra', 'cro', 'pas']

# load pre-processed WETCHIMP results for specified model
# preindustrial partition of wetlands
mod_LSNKcover_save = mod_LSNKcover
for VAR in ['AREA']:
    exec(VAR + '_wet = np.zeros([nb_regionI,nb_biome], dtype=dty)')
    if (mod_EWETpreind != ''):
        for mod_LSNKcover in ['ESA-CCI', 'MODIS', 'Ramankutty1999', 'Levavasseur2012',
                              'mean-TRENDYv2', 'CLM-45', 'JSBACH', 'JULES', 'LPJ',
                              'LPJ-GUESS', 'LPX-Bern', 'OCN', 'ORCHIDEE', 'VISIT']:
            TMP = np.array([line for line in csv.reader(open(
                'data/Wetlands_WETCHIMP/#DATA.Wetlands_' + mod_EWETpreind + '_' + mod_LSNKcover + '.1910s_114reg1_4bio.' + VAR + '.csv',
                'r'))], dtype=dty)
            for i in range(1, 114 + 1):
                for b in range(len(bio) - 2):
                    exec(VAR + '_wet[regionI_index[i],biome_index[bio[b]]] += TMP[i-1,b]')
mod_LSNKcover = mod_LSNKcover_save
del mod_LSNKcover_save

# preindustrial area and emissions
for VAR in ['AREA', 'ECH4']:
    exec(VAR + '_wet0 = np.zeros([nb_regionI], dtype=dty)')
    if (mod_EWETpreind != ''):
        TMP = np.array([line for line in csv.reader(open(
            'data/Wetlands_WETCHIMP/#DATA.Wetlands_' + mod_EWETpreind + '.1910s_114reg1_(1exp).' + VAR + '.csv',
            'r'))][1:], dtype=dty)
        for i in range(1, 114 + 1):
            exec(VAR + '_wet0[regionI_index[i]] += TMP[i-1,0]')

# calculate preindustrial parameters {TgC/Mha}&{Mha}&{.}
# with arbitrary adjustment of areal emissions
ewet_0 = ECH4_wet0 / AREA_wet0 * CO2_0 / np.mean(CO2_cmip5[201:231])
ewet_0[np.isnan(ewet_0) | np.isinf(ewet_0)] = 0
AWET_0 = AREA_wet0.copy()
p_wet = AREA_wet / np.sum(AREA_wet, 1)[:, np.newaxis]
# ensure NaN and zeros removed
for var in ['ewet_0', 'p_wet']:
    exec(var + '[np.isnan(' + var + ')|np.isinf(' + var + ')] = 0')
p_wet[p_wet == 0] = 1E-18
p_wet /= np.sum(p_wet, 1)[:, np.newaxis]

# load pre-processed WETCHIMP results for specified model
# response to perturbation experiments
for sim in ['exp0', 'expC', 'expT', 'expP']:
    for VAR in ['AREA']:
        exec(VAR + '_' + sim + ' = np.zeros([nb_regionI], dtype=dty)')
        if (mod_AWETtrans != ''):
            TMP = np.array([line for line in csv.reader(open(
                'data/Wetlands_WETCHIMP/#DATA.Wetlands_' + mod_AWETtrans + '.1910s_114reg1_(4exp).' + VAR + '.csv',
                'r'))][1:], dtype=dty)
            lgd = [line for line in csv.reader(open(
                'data/Wetlands_WETCHIMP/#DATA.Wetlands_' + mod_AWETtrans + '.1910s_114reg1_(4exp).' + VAR + '.csv',
                'r'))][0]
            for i in range(1, 114 + 1):
                exec(VAR + '_' + sim + '[regionI_index[i]] += TMP[i-1,lgd.index(sim)]')

# value of perturbations
# see [Melton et al., 2013]
CO2_expC = 857. - 303.
T_expT = 3.4
P_expP = np.zeros([nb_regionI], dtype=dty)
S_tmp = np.zeros([nb_regionI], dtype=dty)
TMP = np.array([line for line in csv.reader(
    open('data/HistLand_CRU/#DATA.HistLand_CRU.1901-2014_114reg1.lyp.csv', 'r'))],
               dtype=dty)
TMP2 = np.array([line for line in csv.reader(
    open('data/HistLand_CRU/#DATA.HistLand_CRU.1901-2014_114reg1.AREA.csv', 'r'))],
                dtype=dty)
for i in range(1, 114 + 1):
    P_expP[regionI_index[i]] += np.mean(TMP[:30, i - 1] * TMP2[:30, i - 1], 0)
    S_tmp[regionI_index[i]] += np.mean(TMP2[:30, i - 1], 0)
P_expP /= S_tmp
P_expP[np.isnan(P_expP) | np.isinf(P_expP)] = 0
P_expP *= 0.039

# calculation of parameters
# sensitivity of wetland extent to CO2 {./ppm}
gamma_wetC = (AREA_expC / AREA_exp0 - 1) / CO2_expC
# sensitivity of wetland extent to temperature {./K}
gamma_wetT = (AREA_expT / AREA_exp0 - 1) / T_expT
# sensitivity of wetland extent to precipitation {./mm}
gamma_wetP = (AREA_expP / AREA_exp0 - 1) / P_expP
# [NaN]
for var in ['wetT', 'wetP', 'wetC']:
    exec('gamma_' + var + '[np.isnan(gamma_' + var + ')|np.isinf(gamma_' + var + ')] = 0')

# -----------------
# 2.3.2. Permafrost
# -----------------

# fraction of methane in for PF emissions {.}
# best guess from [Schuur et al., 2015]
if (mod_EPFmethane == 'zero'):
    p_PF_CH4 = np.array([0.000], dtype=dty)
elif (mod_EPFmethane == 'best'):
    p_PF_CH4 = np.array([0.023], dtype=dty)
elif (mod_EPFmethane == 'twice'):
    p_PF_CH4 = np.array([0.046], dtype=dty)

# fraction of instantaneous PF emission {.}
p_PF_inst = np.array([0.0], dtype=dty)

##################################################
#   3. NITROUS OXIDE
##################################################

# ===============
# 3.1. ATMOSPHERE
# ===============

# conversion of N2O from {ppb} to {TgN}
alpha_N2O = 0.1765 * np.array([28.0], dtype=dty)

# historic N2O from IPCC-AR5 {ppb}
# from [IPCC WG1, 2013] annexe 2
N2O_ipcc = np.ones([311 + 1], dtype=dty) * np.nan
TMP = np.array([line for line in csv.reader(
    open('data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1750-2011.N2O.csv', 'r'))],
               dtype=dty)
N2O_ipcc[50:] = TMP[:, 0]
N2O_0 = np.array([N2O_ipcc[50]], dtype=dty)

# historic N2O from CMIP5 {ppb}
# from [Meinshausen et al., 2011]
N2O_cmip5 = np.ones([305 + 1], dtype=dty) * np.nan
TMP = np.array([line for line in csv.reader(
    open('data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005.N2O.csv', 'r'))], dtype=dty)
N2O_cmip5[65:] = TMP[:, 0]

# historic N2O from AGAGE {ppb}
# from [Prinn et al., 2013] updated from the website
N2O_agage = np.ones([313 + 1], dtype=dty) * np.nan
TMP = np.array([line for line in csv.reader(
    open('data/HistAtmo_AGAGE/#DATA.HistAtmo_AGAGE.1979-2013.N2O_global.csv', 'r'))],
               dtype=dty)
N2O_agage[279:] = TMP[:, 0]

# historic N2O from Law Dome ice cores {ppm}
# from [MacFarling Meure et al., 2006]
N2O_lawdome = np.array([line for line in csv.reader(
    open('data/HistAtmo_NOAA-NCDC/#DATA.HistAtmo_NOAA-NCDC.(IceCores).N2O_lawdome.csv',
         'r'))][1:], dtype=dty)

# load RCP concentrations {ppb}
# from [Meinshausen et al., 2011]
N2O_rcp = np.ones([800 + 1, 6], dtype=dty) * np.nan
n = -1
for rcp in ['rcp26', 'rcp45', 'rcp60', 'rcp85', 'rcp45to26', 'rcp60to45']:
    n += 1
    TMP = np.array([line for line in csv.reader(
        open('data/Scenario_ECP/#DATA.Scenario_ECP.2000-2500.' + rcp + '_N2O.csv', 'r'))],
                   dtype=dty)
    N2O_rcp[300:, n] = TMP[:, 0]

# ==============
# 3.2. CHEMISTRY
# ==============

# ----------------
# 3.2.1. Lifetimes
# ----------------

# preindustrial lifetime of N2O {yr}
# average from [Prather et al., 2012]
if (mod_HVSNKtau == 'Prather2015'):
    tau_N2O_hv = np.array([123.], dtype=dty)
# variations also from [Prather et al., 2015] (table 2)
elif (mod_HVSNKtau == 'GMI'):
    tau_N2O_hv = np.array([137.4], dtype=dty) * 123. / 132.5
elif (mod_HVSNKtau == 'GEOSCCM'):
    tau_N2O_hv = np.array([120.2], dtype=dty) * 123. / 132.5
elif (mod_HVSNKtau == 'G2d-M'):
    tau_N2O_hv = np.array([127.0], dtype=dty) * 123. / 132.5
elif (mod_HVSNKtau == 'G2d'):
    tau_N2O_hv = np.array([129.5], dtype=dty) * 123. / 132.5
elif (mod_HVSNKtau == 'Oslo-c29'):
    tau_N2O_hv = np.array([126.1], dtype=dty) * 123. / 132.5
elif (mod_HVSNKtau == 'Oslo-c36'):
    tau_N2O_hv = np.array([146.7], dtype=dty) * 123. / 132.5
elif (mod_HVSNKtau == 'UCI-c29'):
    tau_N2O_hv = np.array([126.2], dtype=dty) * 123. / 132.5
elif (mod_HVSNKtau == 'UCI-c36'):
    tau_N2O_hv = np.array([146.2], dtype=dty) * 123. / 132.5

# lag used for lagged concentrations {yr}
# adjusted; 3yr is typical WMO value
tau_lag = np.array([3.], dtype=dty)

# ------------------
# 3.2.2. Prather2015
# ------------------

# estimate EESC for use to deduce strato sink sensitivity
# ODS concentration from [IPCC WG1, 2013] (annex 2) in year 2005
EESC_hv = np.zeros([nb_ODS], dtype=dty)
EESC_hv0 = np.zeros([nb_ODS], dtype=dty)
for VAR in ODS:
    TMP = np.array([line for line in csv.reader(
        open('data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1960-2011.' + VAR + '.csv',
             'r'))], dtype=dty)
    EESC_hv[ODS.index(VAR)] = TMP[45, :]
EESC_hv0[ODS.index('CH3Br')] = 5.8
EESC_hv0[ODS.index('CH3Cl')] = 480.
#  fractional releases from [Newman et al., 2007]
for VAR in ['EESC_hv', 'EESC_hv0']:
    exec(
        VAR + ' *= np.array([0.47,0.23,0.29,0.12,0.05,0.56,0.67,0.13,0.08,0.01,0.62,0.62,0.28,0.65,0.60,0.44], dtype=dty)')
    exec(
        VAR + ' *= np.array([3,2,3,2,1,4,3,1,2,1,1+60*1,0+60*2,0+60*1,0+60*2,0+60*1,1], dtype=dty)')
    exec(VAR + ' = np.sum(' + VAR + ')')

# sensitivity of strato sink to N2O, ODSs and age-of-air
# from [Prather et al., 2015] (table 3)
# change in age-of-air based on [Fleming et al., 2011] (figure 12)
if (mod_HVSNKtrans == 'Prather2015'):
    chi_hv_N2O = np.array([0.065], dtype=dty)
    chi_hv_EESC = np.array([0.04 / np.log(EESC_hv / EESC_hv0)], dtype=dty)
    chi_hv_age = np.array([0], dtype=dty)
elif (mod_HVSNKtrans == 'G2d'):
    chi_hv_N2O = np.array([0.018 / np.log(321 / 270.)], dtype=dty)
    chi_hv_EESC = np.array([0.048 / np.log(EESC_hv / EESC_hv0)], dtype=dty)
    chi_hv_age = np.array([0.011 / np.log(4.0 / 4.5)], dtype=dty)
elif (mod_HVSNKtrans == 'Oslo-c29'):
    chi_hv_N2O = np.array([0.010 / np.log(321 / 270.)], dtype=dty)
    chi_hv_EESC = np.array([0.033 / np.log(EESC_hv / EESC_hv0)], dtype=dty)
    chi_hv_age = np.array([0], dtype=dty)
elif (mod_HVSNKtrans == 'UCI-c29'):
    chi_hv_N2O = np.array([0.012 / np.log(321 / 270.)], dtype=dty)
    chi_hv_EESC = np.array([0.029 / np.log(EESC_hv / EESC_hv0)], dtype=dty)
    chi_hv_age = np.array([0], dtype=dty)
# but also from [Prather et al., 2012]
elif (mod_HVSNKtrans == 'Prather2012'):
    chi_hv_N2O = np.array([0.08], dtype=dty)
    chi_hv_EESC = np.array([0], dtype=dty)
    chi_hv_age = np.array([0], dtype=dty)


# expression of hv sink function {.}
# adapted from [Prather et al., 2015]
def f_hv(D_N2O, D_EESC, D_gst):
    D_hv = np.exp(chi_hv_N2O * np.log(1 + D_N2O / N2O_0) + chi_hv_EESC * np.log(
        1 + D_EESC / EESC_0) - chi_hv_age * np.log(1 + gamma_age * D_gst)) - 1
    return np.array(D_hv, dtype=dty)


def df_hv_dN2O(D_N2O, D_EESC, D_gst):
    D_hv = chi_hv_N2O / (N2O_0 + D_N2O) * np.exp(
        chi_hv_N2O * np.log(1 + D_N2O / N2O_0) + chi_hv_EESC * np.log(
            1 + D_EESC / EESC_0) - chi_hv_age * np.log(1 + gamma_age * D_gst))
    return np.array(D_hv, dtype=dty)


def df_hv_dEESC(D_N2O, D_EESC, D_gst):
    D_hv = chi_hv_EESC / (EESC_0 + D_EESC) * np.exp(
        chi_hv_N2O * np.log(1 + D_N2O / N2O_0) + chi_hv_EESC * np.log(
            1 + D_EESC / EESC_0) - chi_hv_age * np.log(1 + gamma_age * D_gst))
    return np.array(D_hv, dtype=dty)


def df_hv_dgst(D_N2O, D_EESC, D_gst):
    D_hv = chi_hv_age * gamma_age / (1 + gamma_age * D_gst) * np.exp(
        chi_hv_N2O * np.log(1 + D_N2O / N2O_0) + chi_hv_EESC * np.log(
            1 + D_EESC / EESC_0) - chi_hv_age * np.log(1 + gamma_age * D_gst))
    return np.array(D_hv, dtype=dty)


def df_hv(D_N2O, D_EESC, D_gst):
    df_1 = df_hv_dN2O(D_N2O, D_EESC, D_gst) * D_N2O
    df_2 = df_hv_dEESC(D_N2O, D_EESC, D_gst) * D_EESC
    df_3 = df_hv_dgst(D_N2O, D_EESC, D_gst) * D_gst
    df_tot = df_1 + df_2 + df_3
    return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot),
            np.nan_to_num(df_3 / df_tot)]


# --------------
# 3.2.3. CCMVAL2
# --------------

# load pre-processed CCMVal2 results for specified model
for var in ['age', 'ta']:
    exec(var + '_atm = np.zeros([2099-1961+1], dtype=dty)')
    TMP = np.array([line for line in csv.reader(open(
        'data/Atmosphere_CCMVal2/#DATA.Atmosphere_' + mod_HVSNKcirc + '.1961-2099_(1lvl).' + var + '.csv',
        'r'))][1:], dtype=dty)
    exec(var + '_atm[:] = TMP[:,0]')

# definition of parameter
# sensitivity of stratospheric lag to temperature {./K}
gamma_age = np.array([0], dtype=dty)

# fit of parameter
ratio = np.mean(age_atm[:10]) / age_atm[:]


def err(var):
    clim = 1 + var[0] * (ta_atm[:] - np.mean(ta_atm[:10]))
    return np.sum((ratio - clim) ** 2)


[gamma_age[0]] = fmin(err, [0], disp=False)

##################################################
#   4. HALOGENATED COMPOUNDS
##################################################

# ===============
# 4.1. ATMOSPHERE
# ===============

# conversion of HaloC from {ppt} to {kt}
alpha_HFC = 0.1765 * np.array(
    [70.0, 52.0, 120.0, 102.0, 84.0, 66.0, 170.0, 152.0, 134.0, 148.1, 252.1], dtype=dty)
alpha_PFC = 0.1765 * np.array(
    [146.1, 71.0, 88.0, 138.0, 188.0, 200.0, 238.0, 288.0, 338.0, 388.1], dtype=dty)
alpha_ODS = 0.1765 * np.array(
    [137.4, 120.9, 187.4, 170.9, 154.5, 153.8, 133.4, 86.5, 117.0, 100.5, 165.4, 209.8,
     148.9, 259.8, 94.9, 50.5], dtype=dty)

# historic HaloC from IPCC-AR5 {ppt}
# from [IPCC WG1, 2013] (annexe 2)
for VAR in ['HFC', 'PFC', 'ODS']:
    exec(VAR + '_ipcc = np.ones([311+1,nb_' + VAR + '], dtype=dty) * np.nan')
for VAR in HFC:
    if os.path.isfile(
            'data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1940-2011.' + VAR + '.csv'):
        TMP = np.array([line for line in csv.reader(open(
            'data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1940-2011.' + VAR + '.csv',
            'r'))], dtype=dty)
        HFC_ipcc[240:, HFC.index(VAR)] = TMP[:, 0]
for VAR in PFC:
    if os.path.isfile(
            'data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1900-2011.' + VAR + '.csv'):
        TMP = np.array([line for line in csv.reader(open(
            'data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1900-2011.' + VAR + '.csv',
            'r'))], dtype=dty)
        PFC_ipcc[200:, PFC.index(VAR)] = TMP[:, 0]
for VAR in ODS:
    if os.path.isfile(
            'data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1960-2011.' + VAR + '.csv'):
        TMP = np.array([line for line in csv.reader(open(
            'data/HistAtmo_IPCC-AR5/#DATA.HistAtmo_IPCC-AR5.1960-2011.' + VAR + '.csv',
            'r'))], dtype=dty)
        ODS_ipcc[260:, ODS.index(VAR)] = TMP[:, 0]

# historic HaloC from CMIP5 {ppt}
# from [Meinshausen et al., 2011]
for VAR in ['HFC', 'PFC', 'ODS']:
    exec(VAR + '_cmip5 = np.ones([305+1,nb_' + VAR + '], dtype=dty) * np.nan')
for VAR in HFC:
    if os.path.isfile(
            'data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005.' + VAR + '.csv'):
        TMP = np.array([line for line in csv.reader(
            open('data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005.' + VAR + '.csv',
                 'r'))], dtype=dty)
        HFC_cmip5[65:, HFC.index(VAR)] = TMP[:, 0]
for VAR in PFC:
    if os.path.isfile(
            'data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005.' + VAR + '.csv'):
        TMP = np.array([line for line in csv.reader(
            open('data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005.' + VAR + '.csv',
                 'r'))], dtype=dty)
        PFC_cmip5[65:, PFC.index(VAR)] = TMP[:, 0]
for VAR in ODS:
    if os.path.isfile(
            'data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005.' + VAR + '.csv'):
        TMP = np.array([line for line in csv.reader(
            open('data/HistAtmo_CMIP5/#DATA.HistAtmo_CMIP5.1765-2005.' + VAR + '.csv',
                 'r'))], dtype=dty)
        ODS_cmip5[65:, ODS.index(VAR)] = TMP[:, 0]

# preindustrial HaloC concentrations {ppt}
# from [IPCC WG1, 2013] (annexe 2) and [Meinshausen et al., 2011]
for VAR in ['HFC', 'PFC', 'ODS']:
    exec(VAR + '_0 = np.zeros([nb_' + VAR + '], dtype=dty)')
PFC_0[PFC.index('CF4')] = 35.
ODS_0[ODS.index('CH3Br')] = 5.8
ODS_0[ODS.index('CH3Cl')] = 480.

# ==============
# 4.2. CHEMISTRY
# ==============

# atmospheric OH lifetime {yr}
# from [WMO, 2011] (table 1-3)
# rescaled to follow the arbitrary rescaling of tau_CH4_OH
tau_HFC_OH = scale_OH * np.array([245, 5.5, 32, 14.3, 55, 1.6, 44.5, 253, 8.2, 9.3, 17.9],
                                 dtype=dty)
tau_PFC_OH = scale_OH * np.array(
    [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf],
    dtype=dty)
tau_ODS_OH = scale_OH * np.array(
    [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, 6.1, 12.8, 10.7, 19.3, np.inf,
     np.inf, np.inf, np.inf, 1.5, 1.9], dtype=dty)

# stratospheric lifetime {yr}
# idem
# rescaled to follow [Prather et al., 2015]
tau_HFC_hv = 1.06 * np.array([2347, 89, 246, 232, 327, 45.4, 310, 5676, 116, 125, 157],
                             dtype=dty)
tau_PFC_hv = 1.06 * np.array(
    [3200, 500, 50000, 10000, 2600, 3200, 2600, 4100, 3100, 3000], dtype=dty)
tau_ODS_hv = 1.06 * np.array(
    [45, 100, 85, 190, 1020, 35, 39, 186, 64.9, 160, np.inf, np.inf, np.inf, np.inf,
     np.inf, np.inf], dtype=dty)

# other sinks lifetime {yr}
# idem
tau_HFC_othr = np.array(
    [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf,
     np.inf], dtype=dty)
tau_PFC_othr = np.array(
    [np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf],
    dtype=dty)
tau_ODS_othr = np.array(
    [np.inf, np.inf, np.inf, np.inf, np.inf, 101, 89, np.inf, np.inf, np.inf, 16, 2.9, 65,
     20, 3, 1.4], dtype=dty)

##################################################
#   5. OZONE
##################################################

# ================
# 5.1. TROPOSPHERE
# ================

# -----------
# 5.1.1. HTAP
# -----------

# read region distribution
TMP = np.array([line for line in csv.reader(
    open('data/RegDiv_HTAP/#DATA.RegDiv_HTAP.114reg1_(4reg0).AREA.csv', 'r'))][1:],
               dtype=dty)
p_reg4 = np.zeros([nb_regionI, 4 + 1], dtype=dty)
for i in range(1, 114 + 1):
    p_reg4[regionI_index[i], :] += TMP[i - 1, :]
p_reg4 /= np.sum(p_reg4, 1)[:, np.newaxis]
p_reg4[np.isnan(p_reg4) | np.isinf(p_reg4)] = 0

# regional weight of ozone precursors {.}
# from HTAP experiments [Fry et al., 2013] (table S5) & [Fiore et al., 2013] (table S1)
if (mod_O3Tregsat == 'mean-HTAP'):
    w_reg_NOX = np.array([-3.31, -1.32, -0.67, -0.91, -0.41], dtype=dty) / np.array(
        [27.5, 8.7, 8.4, 7.0, 3.3], dtype=dty) / -0.2
    w_reg_CO = np.array([-1.52, -0.42, -0.27, -0.52, -0.31], dtype=dty) / np.array(
        [461, 123, 88, 151, 98], dtype=dty) / -0.2
    w_reg_VOC = np.array([-2.53, -0.36, -0.42, -0.36, -2.53], dtype=dty) / np.array(
        [191.8, 67.8, 39.1, 49.8, 35.2], dtype=dty) / -0.2
elif (mod_O3Tregsat == 'CAMCHEM'):
    w_reg_NOX = np.array([-1.91, -0.62, -0.44, -0.55, -0.30], dtype=dty) / np.array(
        [27.6, 8.2, 9.2, 6.8, 3.4], dtype=dty) / -0.2
    w_reg_CO = np.array([-1.49, -0.36, -0.22, -0.58, -0.33], dtype=dty) / np.array(
        [626, 154, 106, 221, 145], dtype=dty) / -0.2
    w_reg_VOC = np.array([-1.25, -0.29, -0.59, -0.21, -0.16], dtype=dty) / np.array(
        [245.2, 92.4, 60.6, 55.5, 36.7], dtype=dty) / -0.2
elif (mod_O3Tregsat == 'FRSGCUCI'):
    w_reg_NOX = np.array([-1.58, -0.50, -0.22, -0.62, -0.24], dtype=dty) / np.array(
        [27.1, 8.8, 8.4, 6.9, 3.0], dtype=dty) / -0.2
    w_reg_CO = np.array([-1.54, -0.47, -0.28, -0.49, -0.30], dtype=dty) / np.array(
        [408, 131, 74, 129, 74], dtype=dty) / -0.2
    w_reg_VOC = np.array([-2.42, -0.60, -0.78, -0.68, -0.36], dtype=dty) / np.array(
        [175.1, 57.3, 33.5, 50.5, 33.8], dtype=dty) / -0.2
elif (mod_O3Tregsat == 'GISS-modelE'):
    w_reg_NOX = np.array([-2.84, -1.14, -0.55, -0.71, -0.44], dtype=dty) / np.array(
        [31.4, 8.9, 8.6, 10.8, 3.1], dtype=dty) / -0.2
    w_reg_CO = np.array([-2.23, -0.61, -0.36, -0.90, -0.36], dtype=dty) / np.array(
        [417, 111, 69, 170, 67], dtype=dty) / -0.2
    w_reg_VOC = np.array([-0.41, -0.08, -0.08, -0.14, -0.11], dtype=dty) / np.array(
        [111.7, 29.0, 26.9, 31.7, 24.1], dtype=dty) / -0.2
elif (mod_O3Tregsat == 'GMI'):
    w_reg_NOX = np.array([-2.61, -0.96, -0.54, -0.75, -0.36], dtype=dty) / np.array(
        [24.4, 8.0, 7.5, 5.8, 3.1], dtype=dty) / -0.2
    w_reg_CO = np.array([-2.19, -0.52, -0.38, -0.87, -0.42], dtype=dty) / np.array(
        [528, 134, 101, 194, 99], dtype=dty) / -0.2
    w_reg_VOC = np.array([-1.08, -0.31, -0.15, -0.38, -0.24], dtype=dty) / np.array(
        [165.7, 60.1, 32.0, 42.0, 31.6], dtype=dty) / -0.2
elif (mod_O3Tregsat == 'INCA'):
    w_reg_NOX = np.array([-3.08, -1.13, -0.56, -0.67, -0.72], dtype=dty) / np.array(
        [26.9, 8.6, 7.3, 7.2, 3.8], dtype=dty) / -0.2
    w_reg_CO = np.array([-0.91, -0.18, -0.17, -0.33, -0.23], dtype=dty) / np.array(
        [376, 74, 67, 127, 108], dtype=dty) / -0.2
    w_reg_VOC = np.array([-1.05, -0.40, -0.35, -0.30, 0.00], dtype=dty) / np.array(
        [279.1, 100.1, 56.7, 77.0, 45.3], dtype=dty) / -0.2
elif (mod_O3Tregsat == 'LLNL-IMPACT'):
    w_reg_NOX = np.array([-1.99, -0.82, -0.36, -0.50, -0.31], dtype=dty) / np.array(
        [30.0, 8.8, 9.7, 7.3, 4.2], dtype=dty) / -0.2
    w_reg_CO = np.array([-1.41, -0.31, -0.24, -0.40, -0.46], dtype=dty) / np.array(
        [518, 130, 111, 153, 124], dtype=dty) / -0.2
    w_reg_VOC = np.array([-0.22, -0.05, -0.10, -0.04, -0.03], dtype=dty) / np.array(
        [143.4, 53.2, 23.3, 36.2, 30.7], dtype=dty) / -0.2
elif (mod_O3Tregsat == 'MOZART-GFDL'):
    w_reg_NOX = np.array([-3.00, -1.18, -0.57, -0.81, -0.44], dtype=dty) / np.array(
        [25.8, 9.4, 8.6, 5.2, 2.6], dtype=dty) / -0.2
    w_reg_CO = np.array([-1.29, -0.34, -0.35, -0.36, -0.24], dtype=dty) / np.array(
        [493, 124, 130, 134, 105], dtype=dty) / -0.2
    w_reg_VOC = np.array([-0.84, -0.19, -0.28, -0.20, -0.17], dtype=dty) / np.array(
        [186.8, 68.7, 32.4, 48.9, 36.8], dtype=dty) / -0.2
elif (mod_O3Tregsat == 'MOZECH'):
    w_reg_NOX = np.array([-2.24, -0.81, -0.58, -0.50, -0.35], dtype=dty) / np.array(
        [26.9, 8.9, 7.4, 7.0, 3.6], dtype=dty) / -0.2
    w_reg_CO = np.array([-1.61, -0.46, -0.31, -0.48, -0.36], dtype=dty) / np.array(
        [470, 107, 85, 154, 124], dtype=dty) / -0.2
    w_reg_VOC = np.array([-1.94, -0.69, -0.55, -0.45, -0.25], dtype=dty) / np.array(
        [250.6, 107.0, 44.7, 62.0, 36.9], dtype=dty) / -0.2
elif (mod_O3Tregsat == 'STOC-HadAM3'):
    w_reg_NOX = np.array([-1.43, -0.59, -0.27, -0.46, -0.11], dtype=dty) / np.array(
        [27.9, 9.0, 8.6, 7.1, 3.2], dtype=dty) / -0.2
    w_reg_CO = np.array([-1.29, -0.39, -0.25, -0.40, -0.25], dtype=dty) / np.array(
        [406, 129, 74, 127, 76], dtype=dty) / -0.2
    w_reg_VOC = np.array([-1.76, -0.40, -0.64, -0.45, -0.27], dtype=dty) / np.array(
        [216.9, 70.9, 52.2, 50.0, 43.8], dtype=dty) / -0.2
elif (mod_O3Tregsat == 'TM5-JRC'):
    w_reg_NOX = np.array([-3.80, -1.32, -0.72, -0.98, -0.78], dtype=dty) / np.array(
        [26.6, 8.7, 8.5, 6.4, 3.0], dtype=dty) / -0.2
    w_reg_CO = np.array([-1.34, -0.54, -0.16, -0.42, -0.22], dtype=dty) / np.array(
        [399, 127, 75, 123, 74], dtype=dty) / -0.2
    w_reg_VOC = np.array([-1.91, -0.56, -0.49, -0.58, -0.28], dtype=dty) / np.array(
        [164.3, 50.3, 34.6, 47.0, 32.4], dtype=dty) / -0.2
elif (mod_O3Tregsat == 'UM-CAM'):
    w_reg_NOX = np.array([-2.53, -0.94, -0.50, -0.68, -0.41], dtype=dty) / np.array(
        [27.9, 8.9, 8.7, 7.0, 3.3], dtype=dty) / -0.2
    w_reg_CO = np.array([-1.42, -0.41, -0.28, -0.44, -0.29], dtype=dty) / np.array(
        [428, 137, 81, 131, 79], dtype=dty) / -0.2
    w_reg_VOC = np.array([-1.91, -0.38, -0.64, -0.58, -0.31], dtype=dty) / np.array(
        [171.1, 56.4, 32.7, 47.4, 34.6], dtype=dty) / -0.2
elif (mod_O3Tregsat == ''):
    w_reg_NOX = np.array([1, 1, 1, 1, 1], dtype=dty)
    w_reg_CO = np.array([1, 1, 1, 1, 1], dtype=dty)
    w_reg_VOC = np.array([1, 1, 1, 1, 1], dtype=dty)
w_reg_NOX /= w_reg_NOX[0]
w_reg_CO /= w_reg_CO[0]
w_reg_VOC /= w_reg_VOC[0]

# -------------
# 5.1.2. ACCMIP
# -------------

# load pre-processed ACCMIP results for specified model
# for sensitivity to emissions
if (mod_O3Temis != 'mean-OxComp'):
    for var in ['O3t', 'CH4', 'ECO', 'EVOC', 'ENOX']:
        TMP = np.array([line for line in csv.reader(open(
            'data/OzoChem_ACCMIP/#DATA.OzoChem_' + mod_O3Temis + '.2000s_(6exp).' + var + '.csv',
            'r'))][1:], dtype=dty)
        lgd = [line for line in csv.reader(open(
            'data/OzoChem_ACCMIP/#DATA.OzoChem_' + mod_O3Temis + '.2000s_(6exp).' + var + '.csv',
            'r'))][0]
        exec(var + '_ozoe = TMP[0,:]')

# setting of parameters
# based on ACCMIP
# to be used with formula by [Ehhalt et al., 2001]
if (mod_O3Temis != 'mean-OxComp'):
    # sensitivity of tropospheric O3 to atmospheric CH4 {DU/.}
    chi_O3t_CH4 = np.array(
        [(O3t_ozoe[2] - O3t_ozoe[1]) / np.log(CH4_ozoe[2] / CH4_ozoe[1])], dtype=dty)
    # sensitivity of tropospheric O3 to ozoesions of ozone precursors {DU/{TgC/yr}}&{DU/{Tg/yr}}&{DU/{TgN/yr}}
    chi_O3t_CO = np.array([(O3t_ozoe[3] - O3t_ozoe[1]) / (ECO_ozoe[3] - ECO_ozoe[1])],
                          dtype=dty)
    chi_O3t_VOC = np.array([(O3t_ozoe[4] - O3t_ozoe[1]) / (EVOC_ozoe[4] - EVOC_ozoe[1])],
                           dtype=dty)
    chi_O3t_NOX = np.array([(O3t_ozoe[5] - O3t_ozoe[1]) / (ENOX_ozoe[5] - ENOX_ozoe[1])],
                           dtype=dty)
    # normalization of the non-linearity
    chi_O3t_CH4 *= (O3t_ozoe[1] - O3t_ozoe[0]) / (
                4 * O3t_ozoe[1] - O3t_ozoe[2] - O3t_ozoe[3] - O3t_ozoe[4] - O3t_ozoe[5])
    chi_O3t_CO *= (O3t_ozoe[1] - O3t_ozoe[0]) / (
                4 * O3t_ozoe[1] - O3t_ozoe[2] - O3t_ozoe[3] - O3t_ozoe[4] - O3t_ozoe[5])
    chi_O3t_VOC *= (O3t_ozoe[1] - O3t_ozoe[0]) / (
                4 * O3t_ozoe[1] - O3t_ozoe[2] - O3t_ozoe[3] - O3t_ozoe[4] - O3t_ozoe[5])
    chi_O3t_NOX *= (O3t_ozoe[1] - O3t_ozoe[0]) / (
                4 * O3t_ozoe[1] - O3t_ozoe[2] - O3t_ozoe[3] - O3t_ozoe[4] - O3t_ozoe[5])

# taken from the TAR [Ehhalt et al., 2001]
elif (mod_O3Temis == 'mean-OxComp'):
    chi_O3t_CH4 = np.array([5.0], dtype=dty)
    chi_O3t_CO = np.array([0.0011], dtype=dty) * 28 / 12.
    chi_O3t_VOC = np.array([0.0033], dtype=dty)
    chi_O3t_NOX = np.array([0.125], dtype=dty)

# load pre-processed ACCMIP results for specified model
# for sensitivity to climate
if (mod_O3Tclim != ''):
    for var in ['O3t', 'tas']:
        TMP = np.array([line for line in csv.reader(open(
            'data/OzoChem_ACCMIP/#DATA.OzoChem_' + mod_O3Tclim + '.2000s_(4exp).' + var + '.csv',
            'r'))][1:], dtype=dty)
        exec(var + '_ozoc = TMP[0,:]')

# definition of parameter
# sensitivity of tropospheric O3 to climate change {DU/K}
Gamma_O3t = np.array([0], dtype=dty)

# fit of parameter
if (mod_O3Tclim != ''):
    diff = O3t_ozoc - O3t_ozoc[0]


    def err(var):
        clim = var[0] * (tas_ozoc[:] - tas_ozoc[0])
        return np.sum((diff - clim) ** 2)


    [Gamma_O3t[0]] = fmin(err, [0], disp=False)

# =================
# 5.2. STRATOSPHERE
# =================

# ---------------
# 5.2.1. Chlorine
# ---------------

# fractional release factors {.}
# fits based on [Newman et al., 2006]
if mod_O3Sfracrel in ['Newman2006']:

    def f_fracrel(tau):
        fracrel = np.zeros([nb_ODS], dtype=dty)
        fracrel[ODS.index(
            'CFC11')] = 6.53976E-02 * tau + 4.18938E-02 * tau ** 2 + -3.68985E-03 * tau ** 3
        fracrel[ODS.index(
            'CFC12')] = 4.06080E-02 * tau + 6.74585E-04 * tau ** 2 + 3.70165E-03 * tau ** 3
        fracrel[ODS.index('CFC113')] = (
                                                   5.36055E-02 * tau + 7.16185E-08 * tau ** 2) * np.exp(
            tau / 4.95237E+00)
        fracrel[ODS.index('CFC114')] = (
                                                   2.01017E-02 * tau + -4.67409E-08 * tau ** 2) * np.exp(
            tau / 4.28272E+00)
        fracrel[ODS.index('CFC115')] = (
                                                   3.55000E-03 * tau + 7.67310E-04 * tau ** 2) * np.exp(
            tau / 4.33197E+00)
        fracrel[ODS.index('CCl4')] = (
                                                 8.50117E-02 * tau + 8.33504E-02 * tau ** 2) * np.exp(
            -tau / 5.12976E+00)
        fracrel[ODS.index('CH3CCl3')] = (
                                                    1.63270E-01 * tau + 1.12119E-01 * tau ** 2) * np.exp(
            -tau / 3.74978E+00)
        fracrel[ODS.index(
            'HCFC22')] = 3.02660E-02 * tau + 2.49148E-04 * tau ** 2 + 1.38022E-03 * tau ** 3
        fracrel[ODS.index('HCFC141b')] = (
                                                     2.89263E-03 * tau + 1.14686E-08 * tau ** 2) * np.exp(
            tau / 1.36374E+00)
        fracrel[ODS.index('HCFC142b')] = (
                                                     1.22909E-05 * tau + 1.06509E-04 * tau ** 2) * np.exp(
            tau / 1.23520E+00)
        fracrel[ODS.index('Halon1211')] = (
                                                      -1.51872E-02 * tau + 1.68718E-01 * tau ** 2) * np.exp(
            -tau / 3.43897E+00)
        fracrel[ODS.index('Halon1202')] = (
                                                      -1.51872E-02 * tau + 1.68718E-01 * tau ** 2) * np.exp(
            -tau / 3.43897E+00)
        fracrel[ODS.index('Halon1301')] = (
                                                      5.46396E-02 * tau + 3.03982E-08 * tau ** 2) * np.exp(
            tau / 5.63826E+00)
        fracrel[ODS.index(
            'Halon2402')] = 2.25980E-01 * tau + 2.23266E-04 * tau ** 2 + -1.13882E-03 * tau ** 3
        fracrel[ODS.index('CH3Br')] = (
                                                  7.60148E-02 * tau + 1.12771E-01 * tau ** 2) * np.exp(
            -tau / 4.09494E+00)
        fracrel[ODS.index(
            'CH3Cl')] = 1.39444E-01 * tau + 1.46219E-04 * tau ** 2 + 8.40557E-04 * tau ** 3
        return fracrel


    def df_fracrel_dtau(tau):
        fracrel = np.zeros([nb_ODS], dtype=dty)
        fracrel[ODS.index(
            'CFC11')] = 6.53976E-02 + 2 * 4.18938E-02 * tau + 3 * -3.68985E-03 * tau ** 2
        fracrel[ODS.index(
            'CFC12')] = 4.06080E-02 + 2 * 6.74585E-04 * tau + 3 * 3.70165E-03 * tau ** 2
        fracrel[ODS.index('CFC113')] = (5.36055E-02 + 2 * 7.16185E-08 * tau + (
                    5.36055E-02 / 4.95237E+00) * tau + (
                                                    7.16185E-08 / 4.95237E+00) * tau ** 2) * np.exp(
            tau / 4.95237E+00)
        fracrel[ODS.index('CFC114')] = (2.01017E-02 + 2 * -4.67409E-08 * tau + (
                    2.01017E-02 / 4.28272E+00) * tau + (
                                                    -4.67409E-08 / 4.28272E+00) * tau ** 2) * np.exp(
            tau / 4.28272E+00)
        fracrel[ODS.index('CFC115')] = (3.55000E-03 + 2 * 7.67310E-04 * tau + (
                    3.55000E-03 / 4.33197E+00) * tau + (
                                                    7.67310E-04 / 4.33197E+00) * tau ** 2) * np.exp(
            tau / 4.33197E+00)
        fracrel[ODS.index('CCl4')] = (8.50117E-02 + 2 * 8.33504E-02 * tau - (
                    8.50117E-02 / 5.12976E+00) * tau - (
                                                  8.33504E-02 / 5.12976E+00) * tau ** 2) * np.exp(
            -tau / 5.12976E+00)
        fracrel[ODS.index('CH3CCl3')] = (1.63270E-01 + 2 * 1.12119E-01 * tau - (
                    1.63270E-01 / 3.74978E+00) * tau - (
                                                     1.12119E-01 / 3.74978E+00) * tau ** 2) * np.exp(
            -tau / 3.74978E+00)
        fracrel[ODS.index(
            'HCFC22')] = 3.02660E-02 + 2 * 2.49148E-04 * tau + 3 * 1.38022E-03 * tau ** 2
        fracrel[ODS.index('HCFC141b')] = (2.89263E-03 + 2 * 1.14686E-08 * tau + (
                    2.89263E-03 / 1.36374E+00) * tau + (
                                                      1.14686E-08 / 1.36374E+00) * tau ** 2) * np.exp(
            tau / 1.36374E+00)
        fracrel[ODS.index('HCFC142b')] = (1.22909E-05 + 2 * 1.06509E-04 * tau + (
                    1.22909E-05 / 1.23520E+00) * tau + (
                                                      1.06509E-04 / 1.23520E+00) * tau ** 2) * np.exp(
            tau / 1.23520E+00)
        fracrel[ODS.index('Halon1211')] = (-1.51872E-02 + 2 * 1.68718E-01 * tau - (
                    -1.51872E-02 / 3.43897E+00) * tau - (
                                                       1.68718E-01 / 3.43897E+00) * tau ** 2) * np.exp(
            -tau / 3.43897E+00)
        fracrel[ODS.index('Halon1202')] = (-1.51872E-02 + 2 * 1.68718E-01 * tau - (
                    -1.51872E-02 / 3.43897E+00) * tau - (
                                                       1.68718E-01 / 3.43897E+00) * tau ** 2) * np.exp(
            -tau / 3.43897E+00)
        fracrel[ODS.index('Halon1301')] = (5.46396E-02 + 2 * 3.03982E-08 * tau + (
                    5.46396E-02 / 5.63826E+00) * tau + (
                                                       3.03982E-08 / 5.63826E+00) * tau ** 2) * np.exp(
            tau / 5.63826E+00)
        fracrel[ODS.index(
            'Halon2402')] = 2.25980E-01 + 2 * 2.23266E-04 * tau + 3 * -1.13882E-03 * tau ** 2
        fracrel[ODS.index('CH3Br')] = (7.60148E-02 + 2 * 1.12771E-01 * tau + (
                    7.60148E-02 / 4.09494E+00) * tau + (
                                                   1.12771E-01 / 4.09494E+00) * tau ** 2) * np.exp(
            -tau / 4.09494E+00)
        fracrel[ODS.index(
            'CH3Cl')] = 1.39444E-01 + 2 * 1.46219E-04 * tau + 3 * 8.40557E-04 * tau ** 2
        return fracrel

    # alternative values (mid-latitudes) from [Laube et al., 2013]
# missing values still from [Newman et al., 2006]
elif mod_O3Sfracrel in ['Laube2013']:

    def f_fracrel(tau):
        fracrel = np.zeros([nb_ODS], dtype=dty)
        fracrel[ODS.index('CFC11')] = -0.0173 + 0.098666 * tau + 0.00816955 * tau ** 2
        fracrel[ODS.index('CFC12')] = -0.0154 + 0.046244 * tau + 0.00707356 * tau ** 2
        fracrel[ODS.index('CFC113')] = -0.0059 + 0.049669 * tau + 0.00862413 * tau ** 2
        fracrel[ODS.index('CFC114')] = (
                                                   2.01017E-02 * tau + -4.67409E-08 * tau ** 2) * np.exp(
            tau / 4.28272E+00)  # not given
        fracrel[ODS.index('CFC115')] = (
                                                   3.55000E-03 * tau + 7.67310E-04 * tau ** 2) * np.exp(
            tau / 4.33197E+00)  # not given
        fracrel[ODS.index('CCl4')] = -0.0139 + 0.131338 * tau + 0.00464806 * tau ** 2
        fracrel[ODS.index('CH3CCl3')] = -0.0227 + 0.254820 * tau + -0.01505946 * tau ** 2
        fracrel[ODS.index('HCFC22')] = -0.0190 + 0.021203 * tau + 0.00229434 * tau ** 2
        fracrel[ODS.index('HCFC141b')] = -0.0635 + 0.050362 * tau + 0.00939938 * tau ** 2
        fracrel[ODS.index('HCFC142b')] = -0.0032 + 0.010130 * tau + 0.00233185 * tau ** 2
        fracrel[
            ODS.index('Halon1211')] = -0.0535 + 0.204371 * tau + -0.00464644 * tau ** 2
        fracrel[ODS.index('Halon1202')] = (
                                                      -1.51872E-02 * tau + 1.68718E-01 * tau ** 2) * np.exp(
            -tau / 3.43897E+00)  # not given
        fracrel[ODS.index('Halon1301')] = -0.0185 + 0.061608 * tau + 0.01051828 * tau ** 2
        fracrel[ODS.index(
            'Halon2402')] = 2.25980E-01 * tau + 2.23266E-04 * tau ** 2 + -1.13882E-03 * tau ** 3  # not given
        fracrel[ODS.index('CH3Br')] = (
                                                  7.60148E-02 * tau + 1.12771E-01 * tau ** 2) * np.exp(
            -tau / 4.09494E+00)  # not given
        fracrel[ODS.index(
            'CH3Cl')] = 1.39444E-01 * tau + 1.46219E-04 * tau ** 2 + 8.40557E-04 * tau ** 3  # not given
        return fracrel


    def df_fracrel_dtau(tau):
        fracrel = np.zeros([nb_ODS], dtype=dty)
        fracrel[ODS.index('CFC11')] = 0.098666 + 2 * 0.00816955 * tau
        fracrel[ODS.index('CFC12')] = 0.046244 + 2 * 0.00707356 * tau
        fracrel[ODS.index('CFC113')] = 0.049669 + 2 * 0.00862413 * tau
        fracrel[ODS.index('CFC114')] = (2.01017E-02 + 2 * -4.67409E-08 * tau + (
                    2.01017E-02 / 4.28272E+00) * tau + (
                                                    -4.67409E-08 / 4.28272E+00) * tau ** 2) * np.exp(
            tau / 4.28272E+00)  # not given
        fracrel[ODS.index('CFC115')] = (3.55000E-03 + 2 * 7.67310E-04 * tau + (
                    3.55000E-03 / 4.33197E+00) * tau + (
                                                    7.67310E-04 / 4.33197E+00) * tau ** 2) * np.exp(
            tau / 4.33197E+00)  # not given
        fracrel[ODS.index('CCl4')] = 0.131338 + 2 * 0.00464806 * tau
        fracrel[ODS.index('CH3CCl3')] = 0.254820 + 2 * -0.01505946 * tau
        fracrel[ODS.index('HCFC22')] = 0.021203 + 2 * 0.00229434 * tau
        fracrel[ODS.index('HCFC141b')] = 0.050362 + 2 * 0.00939938 * tau
        fracrel[ODS.index('HCFC142b')] = 0.010130 + 2 * 0.00233185 * tau
        fracrel[ODS.index('Halon1211')] = 0.204371 + 2 * -0.00464644 * tau
        fracrel[ODS.index('Halon1202')] = (-1.51872E-02 + 2 * 1.68718E-01 * tau - (
                    -1.51872E-02 / 3.43897E+00) * tau - (
                                                       1.68718E-01 / 3.43897E+00) * tau ** 2) * np.exp(
            -tau / 3.43897E+00)  # not given
        fracrel[ODS.index('Halon1301')] = 0.061608 + 2 * 0.01051828 * tau
        fracrel[ODS.index(
            'Halon2402')] = 2.25980E-01 + 2 * 2.23266E-04 * tau + 3 * -1.13882E-03 * tau ** 2  # not given
        fracrel[ODS.index('CH3Br')] = (7.60148E-02 + 2 * 1.12771E-01 * tau + (
                    7.60148E-02 / 4.09494E+00) * tau + (
                                                   1.12771E-01 / 4.09494E+00) * tau ** 2) * np.exp(
            -tau / 4.09494E+00)  # not given
        fracrel[ODS.index(
            'CH3Cl')] = 1.39444E-01 + 2 * 1.46219E-04 * tau + 3 * 8.40557E-04 * tau ** 2  # not given
        return fracrel

    # others values for EESC {.}
# relative strength of bromine from [Daniel et al., 2007]
alpha_Br = np.array([60.], dtype=dty)
n_Cl = np.array([3, 2, 3, 2, 1, 4, 3, 1, 2, 1, 1, 0, 0, 0, 0, 1], dtype=dty)
n_Br = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 1, 2, 1, 0], dtype=dty)
EESC_0 = np.sum(f_fracrel(tau_lag) * (n_Cl + alpha_Br * n_Br) * ODS_0)

# --------------
# 5.2.2. CCMVal2
# --------------

# surface temperature trend from CCMVal2
# load pre-processed CCMVal2 results for specified model
for var in ['ta2']:
    exec(var + '_atm = np.zeros([2099-1961+1], dtype=dty)')
    TMP = np.array([line for line in csv.reader(open(
        'data/Atmosphere_CCMVal2/#DATA.Atmosphere_' + mod_O3Strans + '.1961-2099_(1lvl).' + var + '.csv',
        'r'))][1:], dtype=dty)
    exec(var + '_atm[:] = TMP[:,0]')
# fit of linear yearly trend
ta_trend = np.array([0], dtype=dty)
diff = ta2_atm[:] - np.mean(ta2_atm[:10])


def err(var):
    return np.sum((diff - var[0] * np.arange(len(ta2_atm))) ** 2)


[ta_trend[0]] = fmin(err, [0], disp=False)

# sensitivity of strato ozone to chlorine and climate {DU/ppt}&{DU/K}
# from [Douglass et al., 2014] (figure 2; data provided by author)
if (mod_O3Strans == 'mean-CCMVal2'):
    chi_O3s_EESC = np.array([-12.5E-3], dtype=dty)
    Gamma_O3s = np.array([0.012], dtype=dty) / ta_trend
elif (mod_O3Strans == 'AMTRAC'):
    chi_O3s_EESC = np.array([-13.7E-3], dtype=dty)
    Gamma_O3s = np.array([-0.015], dtype=dty) / ta_trend
elif (mod_O3Strans == 'CCSR-NIES'):
    chi_O3s_EESC = np.array([-7.9E-3], dtype=dty)
    Gamma_O3s = np.array([0.013], dtype=dty) / ta_trend
elif (mod_O3Strans == 'CMAM'):
    chi_O3s_EESC = np.array([-8.1E-3], dtype=dty)
    Gamma_O3s = np.array([-0.15], dtype=dty) / ta_trend
elif (mod_O3Strans == 'CNRM-ACM'):
    chi_O3s_EESC = np.array([-23.2E-3], dtype=dty)
    Gamma_O3s = np.array([0.0066], dtype=dty) / ta_trend
elif (mod_O3Strans == 'LMDZrepro'):
    chi_O3s_EESC = np.array([-10.6E-3], dtype=dty)
    Gamma_O3s = np.array([0.042], dtype=dty) / ta_trend
elif (mod_O3Strans == 'MRI'):
    chi_O3s_EESC = np.array([-20.3E-3], dtype=dty)
    Gamma_O3s = np.array([-0.030], dtype=dty) / ta_trend
elif (mod_O3Strans == 'Niwa-SOCOL'):
    chi_O3s_EESC = np.array([-5.2E-3], dtype=dty)
    Gamma_O3s = np.array([0.087], dtype=dty) / ta_trend
elif (mod_O3Strans == 'SOCOL'):
    chi_O3s_EESC = np.array([-10.5E-3], dtype=dty)
    Gamma_O3s = np.array([0.067], dtype=dty) / ta_trend
elif (mod_O3Strans == 'ULAQ'):
    chi_O3s_EESC = np.array([-15.4E-3], dtype=dty)
    Gamma_O3s = np.array([0.062], dtype=dty) / ta_trend
elif (mod_O3Strans == 'UMSLIMCAT'):
    chi_O3s_EESC = np.array([-10.6E-3], dtype=dty)
    Gamma_O3s = np.array([0.0071], dtype=dty) / ta_trend
elif (mod_O3Strans == 'UMUKCA-UCAM'):
    chi_O3s_EESC = np.array([-11.5E-3], dtype=dty)
    Gamma_O3s = np.array([0.04], dtype=dty) / ta_trend

# sensitivity of stratospheric O3 to N2O {DU/ppb}&{DU/ppb/ppt}
# formulation and values from [Daniel et al., 2010]
if (mod_O3Snitrous == 'Daniel2010'):
    chi_O3s_N2O = chi_O3s_EESC * f_fracrel(3)[0] * 6.4 * (1.53 + 0.53 * 240 / 1400.)
    EESC_x = np.array([1400 / 0.53], dtype=dty)
else:
    chi_O3s_N2O = 0 * chi_O3s_EESC
    EESC_x = np.array([np.inf], dtype=dty)

##################################################
#   6. AEROSOLS
##################################################

# ===============
# 6.1. ATMOSPHERE
# ===============

# conversion of SO4 from {TgS} to {Tg(SO4)}
alpha_SO4 = np.array([96 / 32.], dtype=dty)
# conversion of POM from {Tg(OC)} to {Tg(OM)}
if (mod_POAconv == 'default'):
    alpha_POM = np.array([1.4], dtype=dty)
elif (mod_POAconv == 'GFDL'):
    alpha_POM = np.array([1.6], dtype=dty)
elif (mod_POAconv == 'CSIRO'):
    alpha_POM = np.array([1.3], dtype=dty)
# conversion of NO3 from {TgN} to {Tg(NO3)}
alpha_NO3 = np.array([62 / 14.], dtype=dty)

# ==============
# 6.2. CHEMISTRY
# ==============

# -----------
# 6.2.1. HTAP
# -----------

# regional saturation coefficient {.}
# from HTAP experiments [Yu et al., 2013] (table 6 extended; provided by author)
# sulfate
if (mod_SO4regsat == 'mean-HTAP'):
    w_reg_SO2 = np.array([-3.51, -3.87, -3.92, -2.85, -3.92], dtype=dty) / -3.51
elif (mod_SO4regsat == 'CAMCHEM'):
    w_reg_SO2 = np.array([-4.26, -5.01, -4.86, -3.29, -4.55], dtype=dty) / -4.26
elif (mod_SO4regsat == 'GISS-PUCCINI'):
    w_reg_SO2 = np.array([-2.73, -3.88, -2.92, -1.78, -3.44], dtype=dty) / -2.73
elif (mod_SO4regsat == 'GMI'):
    w_reg_SO2 = np.array([-3.61, -3.55, -3.96, -3.30, -3.59], dtype=dty) / -3.61
elif (mod_SO4regsat == 'GOCART'):
    w_reg_SO2 = np.array([-4.05, -3.86, -4.68, -3.79, -3.06], dtype=dty) / -4.05
elif (mod_SO4regsat == 'INCA2'):
    w_reg_SO2 = np.array([-4.03, -4.52, -4.27, -3.33, -5.45], dtype=dty) / -4.03
elif (mod_SO4regsat == 'LLNL-IMPACT'):
    w_reg_SO2 = np.array([-3.70, -3.74, -4.27, -2.84, -4.82], dtype=dty) / -3.70
elif (mod_SO4regsat == 'SPRINTARS'):
    w_reg_SO2 = np.array([-2.16, -2.51, -2.53, -1.64, -2.56], dtype=dty) / -2.16
else:
    w_reg_SO2 = np.array([1, 1, 1, 1, 1], dtype=dty)

# primary organic aerosols
if (mod_POAregsat == 'mean-HTAP'):
    w_reg_OC = np.array([-4.00, -4.39, -4.30, -3.68, -4.09], dtype=dty) / -4.00
elif (mod_POAregsat == 'CAMCHEM'):
    w_reg_OC = np.array([-3.30, -3.90, -3.23, -2.80, -3.71], dtype=dty) / -3.30
elif (mod_POAregsat == 'GISS-PUCCINI'):
    w_reg_OC = np.array([-6.26, -5.86, -6.07, -6.44, -6.38], dtype=dty) / -6.26
elif (mod_POAregsat == 'GMI'):
    w_reg_OC = np.array([-4.62, -5.30, -6.13, -4.22, -4.16], dtype=dty) / -4.62
elif (mod_POAregsat == 'GOCART'):
    w_reg_OC = np.array([-3.57, -3.65, -3.39, -3.28, -4.00], dtype=dty) / -3.57
elif (mod_POAregsat == 'INCA2'):
    w_reg_OC = np.array([-3.07, -4.27, -3.33, -2.88, -2.66], dtype=dty) / -3.07
elif (mod_POAregsat == 'LLNL-IMPACT'):
    w_reg_OC = np.array([-1.31, -1.41, -1.97, -0.99, -1.29], dtype=dty) / -1.31
elif (mod_POAregsat == 'SPRINTARS'):
    w_reg_OC = np.array([-5.86, -6.32, -5.97, -5.12, -6.45], dtype=dty) / -5.86
else:
    w_reg_OC = np.array([1, 1, 1, 1, 1], dtype=dty)

# black carbon
if (mod_BCregsat == 'mean-HTAP'):
    w_reg_BC = np.array([29.51, 27.31, 37.36, 28.36, 25.31], dtype=dty) / 29.51
elif (mod_BCregsat == 'CAMCHEM'):
    w_reg_BC = np.array([27.56, 28.00, 35.71, 25.24, 24.08], dtype=dty) / 27.56
elif (mod_BCregsat == 'GISS-PUCCINI'):
    w_reg_BC = np.array([60.41, 51.67, 69.53, 65.06, 45.46], dtype=dty) / 60.41
elif (mod_BCregsat == 'GMI'):
    w_reg_BC = np.array([26.68, 25.80, 42.13, 24.81, 15.00], dtype=dty) / 26.68
elif (mod_BCregsat == 'GOCART'):
    w_reg_BC = np.array([46.20, 42.30, 52.21, 45.87, 43.71], dtype=dty) / 46.20
elif (mod_BCregsat == 'INCA2'):
    w_reg_BC = np.array([17.32, 16.88, 20.37, 14.14, 23.69], dtype=dty) / 17.32
elif (mod_BCregsat == 'LLNL-IMPACT'):
    w_reg_BC = np.array([7.25, 6.63, 12.99, 5.66, 5.56], dtype=dty) / 7.25
elif (mod_BCregsat == 'SPRINTARS'):
    w_reg_BC = np.array([21.16, 19.93, 28.58, 17.75, 19.69], dtype=dty) / 21.16
else:
    w_reg_BC = np.array([1, 1, 1, 1, 1], dtype=dty)

# -------------
# 6.2.2. ACCMIP
# -------------

# period of ACCMIP/CMIP5 simulations per simulation
prd = {'ctrl': '251yr', 'hist': '1850-2005', 'rcp26': '2006-2100', 'rcp45': '2006-2100',
       'rcp60': '2006-2100', 'rcp85': '2006-2100'}

# load pre-processed ACCMIP/CMIP5 results for specified model
# sulfate
for sim in ['hist', 'rcp26', 'rcp45', 'rcp60', 'rcp85']:
    for VAR in ['ESO2', 'EDMS', 'SO4'] + ['tas']:
        TMP = np.array([line for line in csv.reader(open(
            'data/AeroChem_ACCMIP/#DATA.AeroChem_' + mod_SO4load + '.' + prd[
                sim] + '.' + sim + '_' + VAR + '.csv', 'r'))], dtype=dty)
        exec(VAR + '_' + sim + 'S = TMP[:,0].copy()')
for VAR in ['ESO2', 'EDMS', 'SO4'] + ['tas']:
    exec(
        VAR + '_allS = np.array(list(' + VAR + '_histS)+list(' + VAR + '_histS)+list(' + VAR + '_histS)+list(' + VAR + '_histS)+list(' + VAR + '_rcp26S)+list(' + VAR + '_rcp45S)+list(' + VAR + '_rcp60S)+list(' + VAR + '_rcp85S), dtype=dty)')

# primary organic aerosols
for sim in ['hist', 'rcp26', 'rcp45', 'rcp60', 'rcp85']:
    for VAR in ['EOM', 'EOMBB', 'POA'] + ['tas']:
        TMP = np.array([line for line in csv.reader(open(
            'data/AeroChem_ACCMIP/#DATA.AeroChem_' + mod_POAload + '.' + prd[
                sim] + '.' + sim + '_' + VAR + '.csv', 'r'))], dtype=dty)
        exec(VAR + '_' + sim + 'P = TMP[:,0].copy()')
for VAR in ['EOM', 'EOMBB', 'POA'] + ['tas']:
    exec(
        VAR + '_allP = np.array(list(' + VAR + '_histP)+list(' + VAR + '_histP)+list(' + VAR + '_histP)+list(' + VAR + '_histP)+list(' + VAR + '_rcp26P)+list(' + VAR + '_rcp45P)+list(' + VAR + '_rcp60P)+list(' + VAR + '_rcp85P), dtype=dty)')

# black carbon
for sim in ['hist', 'rcp26', 'rcp45', 'rcp60', 'rcp85']:
    for VAR in ['EBC', 'EBCBB', 'BC'] + ['tas']:
        TMP = np.array([line for line in csv.reader(open(
            'data/AeroChem_ACCMIP/#DATA.AeroChem_' + mod_BCload + '.' + prd[
                sim] + '.' + sim + '_' + VAR + '.csv', 'r'))], dtype=dty)
        exec(VAR + '_' + sim + 'B = TMP[:,0].copy()')
for VAR in ['EBC', 'EBCBB', 'BC'] + ['tas']:
    exec(
        VAR + '_allB = np.array(list(' + VAR + '_histB)+list(' + VAR + '_histB)+list(' + VAR + '_histB)+list(' + VAR + '_histB)+list(' + VAR + '_rcp26B)+list(' + VAR + '_rcp45B)+list(' + VAR + '_rcp60B)+list(' + VAR + '_rcp85B), dtype=dty)')

# nitrate
if not mod_NO3load in ['Bellouin2011', 'Hauglustaine2014']:
    for sim in ['hist', 'rcp26', 'rcp45', 'rcp60', 'rcp85']:
        for VAR in ['ENOX', 'ENH3', 'NO3'] + ['tas2']:
            TMP = np.array([line for line in csv.reader(open(
                'data/AeroChem_ACCMIP/#DATA.AeroChem_' + mod_NO3load + '.' + prd[
                    sim] + '.' + sim + '_' + VAR + '.csv', 'r'))], dtype=dty)
            exec(VAR + '_' + sim + 'N = TMP[:,0].copy()')
    for VAR in ['ENOX', 'ENH3', 'NO3'] + ['tas2']:
        exec(
            VAR + '_allN = np.array(list(' + VAR + '_histN)+list(' + VAR + '_histN)+list(' + VAR + '_histN)+list(' + VAR + '_histN)+list(' + VAR + '_rcp26N)+list(' + VAR + '_rcp45N)+list(' + VAR + '_rcp60N)+list(' + VAR + '_rcp85N), dtype=dty)')

# secondary organic aerosols
if (mod_SOAload != ''):
    for sim in ['hist', 'rcp26', 'rcp45', 'rcp60', 'rcp85']:
        for VAR in ['EVOC', 'EBVOC', 'SOA'] + ['tas2']:
            TMP = np.array([line for line in csv.reader(open(
                'data/AeroChem_ACCMIP/#DATA.AeroChem_' + mod_SOAload + '.' + prd[
                    sim] + '.' + sim + '_' + VAR + '.csv', 'r'))], dtype=dty)
            exec(VAR + '_' + sim + 'Q = TMP[:,0].copy()')
    for VAR in ['EVOC', 'EBVOC', 'SOA'] + ['tas2']:
        exec(
            VAR + '_allQ = np.array(list(' + VAR + '_histQ)+list(' + VAR + '_histQ)+list(' + VAR + '_histQ)+list(' + VAR + '_histQ)+list(' + VAR + '_rcp26Q)+list(' + VAR + '_rcp45Q)+list(' + VAR + '_rcp60Q)+list(' + VAR + '_rcp85Q), dtype=dty)')

# mineral dusts
for sim in ['hist', 'rcp26', 'rcp45', 'rcp60', 'rcp85']:
    for VAR in ['EDUST', 'DUST'] + ['tas']:
        TMP = np.array([line for line in csv.reader(open(
            'data/AeroChem_ACCMIP/#DATA.AeroChem_' + mod_DUSTload + '.' + prd[
                sim] + '.' + sim + '_' + VAR + '.csv', 'r'))], dtype=dty)
        exec(VAR + '_' + sim + 'D = TMP[:,0].copy()')
for VAR in ['EDUST', 'DUST'] + ['tas']:
    exec(
        VAR + '_allD = np.array(list(' + VAR + '_histD)+list(' + VAR + '_histD)+list(' + VAR + '_histD)+list(' + VAR + '_histD)+list(' + VAR + '_rcp26D)+list(' + VAR + '_rcp45D)+list(' + VAR + '_rcp60D)+list(' + VAR + '_rcp85D), dtype=dty)')

# sea salts
for sim in ['hist', 'rcp26', 'rcp45', 'rcp60', 'rcp85']:
    for VAR in ['ESALT', 'SALT'] + ['tas3']:
        TMP = np.array([line for line in csv.reader(open(
            'data/AeroChem_ACCMIP/#DATA.AeroChem_' + mod_SALTload + '.' + prd[
                sim] + '.' + sim + '_' + VAR + '.csv', 'r'))], dtype=dty)
        exec(VAR + '_' + sim + 'T = TMP[:,0].copy()')
for VAR in ['ESALT', 'SALT'] + ['tas3']:
    exec(
        VAR + '_allT = np.array(list(' + VAR + '_histT)+list(' + VAR + '_histT)+list(' + VAR + '_histT)+list(' + VAR + '_histT)+list(' + VAR + '_rcp26T)+list(' + VAR + '_rcp45T)+list(' + VAR + '_rcp60T)+list(' + VAR + '_rcp85T), dtype=dty)')

# definition of parameters
# lifetimes of sulfate precursors {yr}
tau_SO2 = np.array([0], dtype=dty)
tau_DMS = np.array([0], dtype=dty)
# sensitivity of sulfate to climate change {Tg/K}
Gamma_SO4 = np.array([0], dtype=dty)
# lifetimes of primary organic aerosols {yr}
tau_OMff = np.array([0], dtype=dty)
tau_OMbb = np.array([0], dtype=dty)
# sensitivity of primary organic aerosols to climate change {Tg/K}
Gamma_POA = np.array([0], dtype=dty)
# lifetimes of black carbon {yr}
tau_BCff = np.array([0], dtype=dty)
tau_BCbb = np.array([0], dtype=dty)
# sensitivity of black carbon to climate change {Tg/K}
Gamma_BC = np.array([0], dtype=dty)
# lifetimes of nitrate precursors {yr}
tau_NOX = np.array([0], dtype=dty)
tau_NH3 = np.array([0], dtype=dty)
# sensitivity of nitrate to climate change {Tg/K}
Gamma_NO3 = np.array([0], dtype=dty)
# lifetimes of secondary organic aerosols {yr}
tau_VOC = np.array([0], dtype=dty)
tau_BVOC = np.array([0], dtype=dty)
# sensitivity of secondary organic aerosols to climate change {Tg/K}
Gamma_SOA = np.array([0], dtype=dty)
# lifetime of mineral dusts {yr}
tau_DUST = np.array([0], dtype=dty)
# sensitivity of mineral dusts to climate change {Tg/K}
Gamma_DUST = np.array([0], dtype=dty)
# lifetime of sea salts {yr}
tau_SALT = np.array([0], dtype=dty)
# sensitivity of sea salts to climate change {Tg/K}
Gamma_SALT = np.array([0], dtype=dty)

# fit of parameters
# sulfate
diff = SO4_allS - np.mean(SO4_allS[:10])


def err(var):
    conc = np.abs(var[0]) * alpha_SO4 * (ESO2_allS - np.mean(ESO2_allS[:10])) + np.abs(
        var[1]) * alpha_SO4 * (EDMS_allS - np.mean(EDMS_allS[:10]))
    clim = var[2] * (tas_allS - np.mean(tas_allS[:10]))
    return np.sum((diff - (conc + clim)) ** 2)


[tau_SO2[0], tau_DMS[0], Gamma_SO4[0]] = fmin(err, [0.01, 0.01, 0], disp=False)
tau_SO2 = np.abs(tau_SO2)
tau_DMS = np.abs(tau_DMS)

# primary organic aerosols
diff = POA_allP - np.mean(POA_allP[:10])


def err(var):
    conc = np.abs(var[0]) * (EOM_allP - np.mean(EOM_allP[:10])) + np.abs(var[1]) * (
                EOMBB_allP - np.mean(EOMBB_allP[:10]))
    clim = var[2] * (tas_allP - np.mean(tas_allP[:10]))
    return np.sum((diff - (conc + clim)) ** 2)


[tau_OMff[0], tau_OMbb[0], Gamma_POA[0]] = fmin(err, [0.01, 0.01, 0], disp=False)
tau_OMff = np.abs(tau_OMff)
tau_OMbb = np.abs(tau_OMbb)

# black carbon
diff = BC_allB - np.mean(BC_allB[:10])


def err(var):
    conc = np.abs(var[0]) * (EBC_allB - np.mean(EBC_allB[:10])) + np.abs(var[1]) * (
                EBCBB_allB - np.mean(EBCBB_allB[:10]))
    clim = var[2] * (tas_allB - np.mean(tas_allB[:10]))
    return np.sum((diff - (conc + clim)) ** 2)


[tau_BCff[0], tau_BCbb[0], Gamma_BC[0]] = fmin(err, [0.01, 0.01, 0], disp=False)
tau_BCff = np.abs(tau_BCff)
tau_BCbb = np.abs(tau_BCbb)

# nitrate
if not mod_NO3load in ['Bellouin2011', 'Hauglustaine2014']:
    diff = NO3_allN - np.mean(NO3_allN[:10])


    def err(var):
        conc = np.abs(var[0]) * alpha_NO3 * (
                    ENOX_allN - np.mean(ENOX_allN[:10])) + np.abs(var[1]) * alpha_NO3 * (
                           ENH3_allN - np.mean(ENH3_allN[:10]))
        clim = var[2] * (tas2_allN - np.mean(tas2_allN[:10]))
        return np.sum((diff - (conc + clim)) ** 2)


    [tau_NOX[0], tau_NH3[0], Gamma_NO3[0]] = fmin(err, [0.01, 0.01, 0], disp=False)
    tau_NOX = np.abs(tau_NOX)
    tau_NH3 = np.abs(tau_NH3)

# secondary organic aerosols
if (mod_SOAload != ''):
    diff = SOA_allQ - np.mean(SOA_allQ[:10])


    def err(var):
        conc = np.abs(var[0]) * (EVOC_allQ - np.mean(EVOC_allQ[:10])) + np.abs(var[1]) * (
                    EBVOC_allQ - np.mean(EBVOC_allQ[:10]))
        clim = var[2] * (tas2_allQ - np.mean(tas2_allQ[:10]))
        return np.sum((diff - (conc + clim)) ** 2)


    [tau_VOC[0], tau_BVOC[0], Gamma_SOA[0]] = fmin(err, [0.01, 0.01, 0], disp=False)
    tau_VOC = np.abs(tau_VOC)
    tau_BVOC = np.abs(tau_BVOC)

# mineral dusts
diff = DUST_allD - np.mean(DUST_allD[:10])


def err(var):
    conc = np.abs(var[0]) * (EDUST_allD - np.mean(EDUST_allD[:10]))
    clim = var[1] * (tas_allD - np.mean(tas_allD[:10]))
    return np.sum((diff - (conc + clim)) ** 2)


[tau_DUST[0], Gamma_DUST[0]] = fmin(err, [0.01, 0], disp=False)
tau_DUST = np.abs(tau_DUST)

# sea salts
diff = SALT_allT - np.mean(SALT_allT[:10])


def err(var):
    conc = np.abs(var[0]) * (ESALT_allT - np.mean(ESALT_allT[:10]))
    clim = var[1] * (tas3_allT - np.mean(tas3_allT[:10]))
    return np.sum((diff - (conc + clim)) ** 2)


[tau_SALT[0], Gamma_SALT[0]] = fmin(err, [0.01, 0], disp=False)
tau_SALT = np.abs(tau_SALT)

# ---------------
# 6.2.3. Nitrates
# ---------------

# load data for nitrate aerosols
# from HadGEM2 [Bellouin et al., 2011] (also RCP and CMIP5 data)
if (mod_NO3load == 'Bellouin2011'):
    ENOX_nitrate = np.array([5.7, 37.4, 18.4, 16.3, 16.6, 23.8, 18.4, 16.3, 16.6, 23.8],
                            dtype=dty)
    ENH3_nitrate = np.array([16.6, 41.4, 67.2, 49.2, 63.0, 70.0, 67.2, 49.2, 63.0, 70.0],
                            dtype=dty)
    tas_nitrate = np.array(
        [13.55, 14.07, 15.39, 16.46, 17.07, 18.66, 13.55, 13.55, 13.55, 13.55], dtype=dty)
    NO3_nitrate = np.array([0.05, 0.34, 0.56, 0.29, 0.41, 0.52, 0.63, 0.36, 0.50, 0.68],
                           dtype=dty)
# from LMDz4-INCA3 [Hauglustaine et al., 2014] (tables 1,5)
elif (mod_NO3load == 'Hauglustaine2014'):
    ENOX_nitrate = np.array(
        [10, 36, 29, 26, 14, 32, 26, 14, 30, 27, 13, 38, 30, 21, 21, 14], dtype=dty)
    ENH3_nitrate = np.array(
        [21, 29, 41, 46, 58, 35, 36, 33, 36, 43, 51, 42, 48, 57, 33, 57], dtype=dty)
    tas_nitrate = np.array(
        [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
         0.00, 0.00, 0.00], dtype=dty)
    # NO3_nitrate = np.array([0.23,0.48,0.46,0.47,0.37,0.48,0.46,0.42,0.47,0.48,0.40,0.54,0.52,0.52,0.47,0.43], dtype=dty) # HNO3 and NO3-
    NO3_nitrate = np.array(
        [0.09, 0.18, 0.21, 0.23, 0.21, 0.2, 0.2, 0.18, 0.19, 0.21, 0.21, 0.23, 0.24, 0.25,
         0.2, 0.22], dtype=dty)  # NO3- only

# fit of parameters (defined in previous section)
diff = NO3_nitrate - NO3_nitrate[0]


def err(var):
    conc = np.abs(var[0]) * alpha_NO3 * (ENOX_nitrate - ENOX_nitrate[0]) + np.abs(
        var[1]) * alpha_NO3 * (ENH3_nitrate - ENH3_nitrate[0])
    clim = var[2] * (tas_nitrate - tas_nitrate[0])
    if (mod_NO3load == 'Hauglustaine2014'):
        clim = 0 * clim
    return np.sum((diff - (conc + clim)) ** 2)


[tau_NOX[0], tau_NH3[0], Gamma_NO3[0]] = fmin(err, [0.01, 0.01, 0], disp=False)
tau_NOX = np.abs(tau_NOX)
tau_NH3 = np.abs(tau_NH3)
if (mod_NO3load == 'Hauglustaine2014'):
    Gamma_NO3 = 0 * Gamma_NO3

##################################################
#   7. RADIATIVE FORCING
##################################################

# ====================
# 7.A. Reconstructions
# ====================

# historic RF from IPCC-AR5 {W/m2}
# from [IPCC WG1, 2013] annexe 2
RF_ipcc = np.zeros([311 + 1], dtype=dty)
RF_ipcc[:50] = np.nan
for VAR in ['WMGHG', 'O3', 'AER', 'Alb']:
    exec('RF_' + VAR + '_ipcc = RF_ipcc.copy()')
TMP = np.array([line for line in csv.reader(
    open('data/Historic_IPCC-AR5/#DATA.Historic_IPCC-AR5.1750-2011_(11for).RF.csv',
         'r'))][1:], dtype=dty)
lgd = [line for line in csv.reader(
    open('data/Historic_IPCC-AR5/#DATA.Historic_IPCC-AR5.1750-2011_(11for).RF.csv',
         'r'))][0]
for x in range(len(lgd)):
    if lgd[x] in ['CO2', 'GHG Other', 'H2O (Strat)']:
        RF_WMGHG_ipcc[51:] += TMP[1:, x]
        RF_ipcc[51:] += TMP[1:, x]
    elif lgd[x] in ['O3 (Trop)', 'O3 (Strat)']:
        RF_O3_ipcc[51:] += TMP[1:, x]
        RF_ipcc[51:] += TMP[1:, x]
    elif lgd[x] in ['Aerosol (Total)']:
        RF_AER_ipcc[51:] += TMP[1:, x]
        RF_ipcc[51:] += TMP[1:, x]
    elif lgd[x] in ['LUC', 'BC Snow']:
        RF_Alb_ipcc[51:] += TMP[1:, x]
        RF_ipcc[51:] += TMP[1:, x]
    elif lgd[x] in ['Volcano']:
        RF_ipcc[51:] += TMP[1:, x] - np.mean(TMP[1:, x])
    else:
        RF_ipcc[51:] += TMP[1:, x]

# historic RF from CMIP5 {W/m2}
# from [Meinshausen et al., 2011]
RF_cmip5 = np.zeros([305 + 1], dtype=dty)
RF_cmip5[:65] = np.nan
TMP = np.array([line for line in csv.reader(
    open('data/Historic_CMIP5/#DATA.Historic_CMIP5.1765-2005_(19for).RF.csv', 'r'))][1:],
               dtype=dty)
lgd = [line for line in csv.reader(
    open('data/Historic_CMIP5/#DATA.Historic_CMIP5.1765-2005_(19for).RF.csv', 'r'))][0]
for x in range(len(lgd)):
    if (lgd[x] == 'VOLC'):
        RF_cmip5[66:] += TMP[1:, x] - np.mean(TMP[1:, x])
    else:
        RF_cmip5[66:] += TMP[1:, x]

# load RCP radiative forcing {W/m2}
# from [Meinshausen et al., 2011]
RF_rcp = np.zeros([800 + 1, 6], dtype=dty)
RF_rcp[:300] = np.nan
n = -1
for rcp in ['rcp26', 'rcp45', 'rcp60', 'rcp85', 'rcp45to26', 'rcp60to45']:
    n += 1
    TMP = np.array([line for line in csv.reader(
        open('data/Scenario_ECP/#DATA.Scenario_ECP.2000-2500_(19for).' + rcp + '_RF.csv',
             'r'))][1:], dtype=dty)
    TMP2 = np.array([line for line in csv.reader(
        open('data/Historic_CMIP5/#DATA.Historic_CMIP5.1765-2005_(19for).RF.csv', 'r'))][
                    1:], dtype=dty)
    lgd = [line for line in csv.reader(
        open('data/Scenario_ECP/#DATA.Scenario_ECP.2000-2500_(19for).' + rcp + '_RF.csv',
             'r'))][0]
    for x in range(len(lgd)):
        RF_rcp[300:, n] += TMP[:, x]
        if (lgd[x] == 'VOLC'):
            RF_rcp[:7] -= np.mean(TMP2[1:, x])


# =====================
# 7.1. GREENHOUSE GASES
# =====================

# radiative forcing functions {W/m2}
# from IPCC-TAR [Ramaswamy et al., 2001]
def f_RF_CO2(D_CO2):
    RF = 5.35 * np.log(1 + D_CO2 / CO2_0)
    return np.array(RF, dtype=dty)


def f_RF_CH4(D_CH4):
    RF = 0.036 * (np.sqrt(CH4_0 + D_CH4) - np.sqrt(CH4_0))
    return np.array(RF, dtype=dty)


def f_RF_H2Os(D_CH4_lag):
    RF = 0.15 * 0.036 * (np.sqrt(CH4_0 + D_CH4_lag) - np.sqrt(CH4_0))
    return np.array(RF, dtype=dty)


def f_RF_N2O(D_N2O):
    RF = 0.12 * (np.sqrt(N2O_0 + D_N2O) - np.sqrt(N2O_0))
    return np.array(RF, dtype=dty)


def f_RF_overlap(D_CH4, D_N2O):
    RF = 0.47 * np.log(
        1 + 2.01E-5 * ((CH4_0 + D_CH4) * (N2O_0 + D_N2O)) ** 0.75 + 5.31E-15 * (
                    CH4_0 + D_CH4) * ((CH4_0 + D_CH4) * (N2O_0 + D_N2O)) ** 1.52)
    RF -= 0.47 * np.log(1 + 2.01E-5 * (CH4_0 * N2O_0) ** 0.75 + 5.31E-15 * CH4_0 * (
                CH4_0 * N2O_0) ** 1.52)
    return np.array(RF, dtype=dty)


def df_RF_overlap_dCH4(D_CH4, D_N2O):
    RF = 0.47 * (2.01E-5 * 0.75 * (CH4_0 + D_CH4) ** (0.75 - 1) * (N2O_0 + D_N2O) ** (
        0.75) + 5.31E-15 * (1.52 + 1) * ((CH4_0 + D_CH4) * (N2O_0 + D_N2O)) ** 1.52)
    RF /= 1 + 2.01E-5 * ((CH4_0 + D_CH4) * (N2O_0 + D_N2O)) ** 0.75 + 5.31E-15 * (
                CH4_0 + D_CH4) * ((CH4_0 + D_CH4) * (N2O_0 + D_N2O)) ** 1.52
    return np.array(RF, dtype=dty)


def df_RF_overlap_dN2O(D_CH4, D_N2O):
    RF = 0.47 * (2.01E-5 * 0.75 * (CH4_0 + D_CH4) ** (0.75) * (N2O_0 + D_N2O) ** (
                0.75 - 1) + 5.31E-15 * 1.52 * (CH4_0 + D_CH4) ** (1.52 + 1) * (
                             N2O_0 + D_N2O) ** (1.52 - 1))
    RF /= 1 + 2.01E-5 * ((CH4_0 + D_CH4) * (N2O_0 + D_N2O)) ** 0.75 + 5.31E-15 * (
                CH4_0 + D_CH4) * ((CH4_0 + D_CH4) * (N2O_0 + D_N2O)) ** 1.52
    return np.array(RF, dtype=dty)


def df_RF_overlap(D_CH4, D_N2O):
    df_1 = df_RF_overlap_dCH4(D_CH4, D_N2O) * D_CH4
    df_2 = df_RF_overlap_dN2O(D_CH4, D_N2O) * D_N2O
    df_tot = df_1 + df_2
    return [np.nan_to_num(df_1 / df_tot), np.nan_to_num(df_2 / df_tot)]


# radiative efficiency of halogenated compounds {{W/m2}/ppt}
# from IPCC-AR5 [Myrhe et al., 2013] (table 8.A.1)
radeff_HFC = 1E-3 * np.array(
    [0.18, 0.11, 0.23, 0.16, 0.16, 0.10, 0.26, 0.24, 0.24, 0.22, 0.42], dtype=dty)
radeff_PFC = 1E-3 * np.array([0.57, 0.20, 0.09, 0.25, 0.28, 0.32, 0.36, 0.41, 0.44, 0.50],
                             dtype=dty)
radeff_ODS = 1E-3 * np.array(
    [0.26, 0.32, 0.30, 0.31, 0.20, 0.17, 0.07, 0.21, 0.16, 0.19, 0.29, 0.27, 0.30, 0.31,
     0.004, 0.01], dtype=dty)

# ==========
# 7.2. OZONE
# ==========

# radiative efficiency of tropospheric O3 {{W/m2}/DU}
# from IPCC-AR5 [Myhre et al., 2013]
if (mod_O3Tradeff == 'IPCC-AR5'):
    radeff_O3t = np.array([0.042], dtype=dty)
# from IPCC-AR4 [Forster et al., 2007]
elif (mod_O3Tradeff == 'IPCC-AR4'):
    radeff_O3t = np.array([0.032], dtype=dty)
# from ACCMIP [Stevenson et al., 2013] (table 3)
elif (mod_O3Tradeff == 'mean-ACCMIP'):
    radeff_O3t = np.array([0.377 / 8.9], dtype=dty)
elif (mod_O3Tradeff == 'CESM-CAM-superfast'):
    radeff_O3t = np.array([0.446 / 10.0], dtype=dty)
elif (mod_O3Tradeff == 'CICERO-OsloCTM2'):
    radeff_O3t = np.array([0.401 / 9.3], dtype=dty)
elif (mod_O3Tradeff == 'CMAM'):
    radeff_O3t = np.array([0.322 / 7.6], dtype=dty)
elif (mod_O3Tradeff == 'EMAC'):
    radeff_O3t = np.array([0.460 / 10.8], dtype=dty)
elif (mod_O3Tradeff == 'GEOSCCM'):
    radeff_O3t = np.array([0.387 / 8.7], dtype=dty)
elif (mod_O3Tradeff == 'GFDL-AM3'):
    radeff_O3t = np.array([0.423 / 10.3], dtype=dty)
elif (mod_O3Tradeff == 'GISS-E2-R'):
    radeff_O3t = np.array([0.314 / 8.3], dtype=dty)
elif (mod_O3Tradeff == 'GISS-E2-R-TOMAS'):
    radeff_O3t = np.array([0.333 / 8.7], dtype=dty)
elif (mod_O3Tradeff == 'HadGEM2'):
    radeff_O3t = np.array([0.303 / 7.3], dtype=dty)
elif (mod_O3Tradeff == 'LMDzORINCA'):
    radeff_O3t = np.array([0.351 / 8.2], dtype=dty)
elif (mod_O3Tradeff == 'MIROC-CHEM'):
    radeff_O3t = np.array([0.402 / 9.2], dtype=dty)
elif (mod_O3Tradeff == 'MOCAGE'):
    radeff_O3t = np.array([0.219 / 4.8], dtype=dty)
elif (mod_O3Tradeff == 'NCAR-CAM-35'):
    radeff_O3t = np.array([0.433 / 10.2], dtype=dty)
elif (mod_O3Tradeff == 'STOC-HadAM3'):
    radeff_O3t = np.array([0.437 / 10.5], dtype=dty)
elif (mod_O3Tradeff == 'UM-CAM'):
    radeff_O3t = np.array([0.376 / 8.7], dtype=dty)
elif (mod_O3Tradeff == 'TM5'):
    radeff_O3t = np.array([0.422 / 10.0], dtype=dty)

# radiative efficiency of stratospheric O3 {{W/m2}/DU}
# from IPCC-AR4 [Forster et al., 2007]
if (mod_O3Sradeff == 'IPCC-AR4'):
    radeff_O3s = np.array([0.004], dtype=dty)
# from ACCENT [Gauss et al., 2006] (tables 4 & 6)
elif (mod_O3Sradeff == 'mean-ACCENT'):
    radeff_O3s = np.array([-0.058 / -13.9], dtype=dty)
elif (mod_O3Sradeff == 'ULAQ'):
    radeff_O3s = np.array([-0.059 / -12.6], dtype=dty)
elif (mod_O3Sradeff == 'DLR-E39C'):
    radeff_O3s = np.array([-0.027 / -16.1], dtype=dty)
elif (mod_O3Sradeff == 'NCAR-MACCM'):
    radeff_O3s = np.array([-0.019 / -12.7], dtype=dty)
elif (mod_O3Sradeff == 'CHASER'):
    radeff_O3s = np.array([-0.126 / -14.1], dtype=dty)

# =============
# 7.3. AEROSOLS
# =============

# -------------
# 7.3.1. Direct
# -------------

# radiative efficiency of sulfate aerosols {{W/m2}/Tg}
# from AeroCom2 [Myhre et al., 2013] (table 4)
if (mod_SO4radeff == 'mean-AeroCom2'):
    radeff_SO4 = 1E12 / 510072E9 * np.array([-185], dtype=dty)
elif (mod_SO4radeff == 'BCC'):
    radeff_SO4 = 1E12 / 510072E9 * np.array([-108], dtype=dty)
elif (mod_SO4radeff == 'CAM4-Oslo'):
    radeff_SO4 = 1E12 / 510072E9 * np.array([-173], dtype=dty)
elif (mod_SO4radeff == 'CAM-51'):
    radeff_SO4 = 1E12 / 510072E9 * np.array([-104], dtype=dty)
elif (mod_SO4radeff == 'GEOS-CHEM'):
    radeff_SO4 = 1E12 / 510072E9 * np.array([-123], dtype=dty)
elif (mod_SO4radeff == 'GISS-MATRIX'):
    radeff_SO4 = 1E12 / 510072E9 * np.array([-196], dtype=dty)
elif (mod_SO4radeff == 'GISS-modelE'):
    radeff_SO4 = 1E12 / 510072E9 * np.array([-307], dtype=dty)
elif (mod_SO4radeff == 'GMI'):
    radeff_SO4 = 1E12 / 510072E9 * np.array([-195], dtype=dty)
elif (mod_SO4radeff == 'GOCART'):
    radeff_SO4 = 1E12 / 510072E9 * np.array([-238], dtype=dty)
elif (mod_SO4radeff == 'HadGEM2'):
    radeff_SO4 = 1E12 / 510072E9 * np.array([-193], dtype=dty)
elif (mod_SO4radeff == 'IMPACT-Umich'):
    radeff_SO4 = 1E12 / 510072E9 * np.array([-113], dtype=dty)
elif (mod_SO4radeff == 'INCA'):
    radeff_SO4 = 1E12 / 510072E9 * np.array([-180], dtype=dty)
elif (mod_SO4radeff == 'MPIHAM'):
    radeff_SO4 = 1E12 / 510072E9 * np.array([-125], dtype=dty)
elif (mod_SO4radeff == 'NCAR-CAM-35'):
    radeff_SO4 = 1E12 / 510072E9 * np.array([-354], dtype=dty)
elif (mod_SO4radeff == 'OsloCTM2'):
    radeff_SO4 = 1E12 / 510072E9 * np.array([-192], dtype=dty)
elif (mod_SO4radeff == 'SPRINTARS'):
    radeff_SO4 = 1E12 / 510072E9 * np.array([-172], dtype=dty)

# radiative efficiency of primary organic aerosols {{W/m2}/Tg}
# from AeroCom2 [Myhre et al., 2013] (table 6)
if (mod_POAradeff == 'mean-AeroCom2'):
    radeff_POA = 1E12 / 510072E9 * np.array([-113], dtype=dty)
elif (mod_POAradeff == 'BCC'):
    radeff_POA = 1E12 / 510072E9 * np.array([-97], dtype=dty)
elif (mod_POAradeff == 'CAM4-Oslo'):
    radeff_POA = 1E12 / 510072E9 * np.array([-118], dtype=dty)
elif (mod_POAradeff == 'CAM-51'):
    radeff_POA = 1E12 / 510072E9 * np.array([-69], dtype=dty)
elif (mod_POAradeff == 'GEOS-CHEM'):
    radeff_POA = 1E12 / 510072E9 * np.array([-95], dtype=dty)
elif (mod_POAradeff == 'GISS-MATRIX'):
    radeff_POA = 1E12 / 510072E9 * np.array([-129], dtype=dty)
elif (mod_POAradeff == 'GISS-modelE'):
    radeff_POA = 1E12 / 510072E9 * np.array([-76], dtype=dty)
elif (mod_POAradeff == 'GMI'):
    radeff_POA = 1E12 / 510072E9 * np.array([-189], dtype=dty)
elif (mod_POAradeff == 'GOCART'):
    radeff_POA = 1E12 / 510072E9 * np.array([-144], dtype=dty)
elif (mod_POAradeff == 'HadGEM2'):
    radeff_POA = 1E12 / 510072E9 * np.array([-145], dtype=dty)
elif (mod_POAradeff == 'IMPACT-Umich'):
    radeff_POA = 1E12 / 510072E9 * np.array([-141], dtype=dty)
elif (mod_POAradeff == 'INCA'):
    radeff_POA = 1E12 / 510072E9 * np.array([-76], dtype=dty)
elif (mod_POAradeff == 'MPIHAM'):
    radeff_POA = 1E12 / 510072E9 * np.array([-41], dtype=dty)
elif (mod_POAradeff == 'NCAR-CAM-35'):
    radeff_POA = 1E12 / 510072E9 * np.array([-48], dtype=dty)
elif (mod_POAradeff == 'OsloCTM2'):
    radeff_POA = 1E12 / 510072E9 * np.array([-165], dtype=dty)
elif (mod_POAradeff == 'SPRINTARS'):
    radeff_POA = 1E12 / 510072E9 * np.array([-102], dtype=dty)

# radiative efficiency of black carbon aerosols {{W/m2}/Tg}
# from AeroCom2 [Myhre et al., 2013] (table 5)
if (mod_BCradeff == 'mean-AeroCom2'):
    radeff_BC = 1E12 / 510072E9 * np.array([1438], dtype=dty)
elif (mod_BCradeff == 'BCC'):
    radeff_BC = 1E12 / 510072E9 * np.array([650], dtype=dty)
elif (mod_BCradeff == 'CAM4-Oslo'):
    radeff_BC = 1E12 / 510072E9 * np.array([1763], dtype=dty)
elif (mod_BCradeff == 'CAM-51'):
    radeff_BC = 1E12 / 510072E9 * np.array([2661], dtype=dty)
elif (mod_BCradeff == 'GEOS-CHEM'):
    radeff_BC = 1E12 / 510072E9 * np.array([1067], dtype=dty)
elif (mod_BCradeff == 'GISS-MATRIX'):
    radeff_BC = 1E12 / 510072E9 * np.array([2484], dtype=dty)
elif (mod_BCradeff == 'GISS-modelE'):
    radeff_BC = 1E12 / 510072E9 * np.array([1253], dtype=dty)
elif (mod_BCradeff == 'GMI'):
    radeff_BC = 1E12 / 510072E9 * np.array([1208], dtype=dty)
elif (mod_BCradeff == 'GOCART'):
    radeff_BC = 1E12 / 510072E9 * np.array([874], dtype=dty)
elif (mod_BCradeff == 'HadGEM2'):
    radeff_BC = 1E12 / 510072E9 * np.array([612], dtype=dty)
elif (mod_BCradeff == 'IMPACT-Umich'):
    radeff_BC = 1E12 / 510072E9 * np.array([1467], dtype=dty)
elif (mod_BCradeff == 'INCA'):
    radeff_BC = 1E12 / 510072E9 * np.array([1160], dtype=dty)
elif (mod_BCradeff == 'MPIHAM'):
    radeff_BC = 1E12 / 510072E9 * np.array([1453], dtype=dty)
elif (mod_BCradeff == 'NCAR-CAM-35'):
    radeff_BC = 1E12 / 510072E9 * np.array([1364], dtype=dty)
elif (mod_BCradeff == 'OsloCTM2'):
    radeff_BC = 1E12 / 510072E9 * np.array([2161], dtype=dty)
elif (mod_BCradeff == 'SPRINTARS'):
    radeff_BC = 1E12 / 510072E9 * np.array([1322], dtype=dty)

# radiative efficiency of nitrate aerosols {{W/m2}/Tg}
# from AeroCom2 [Myhre et al., 2013] (table 8)
if (mod_NO3radeff == 'mean-AeroCom2'):
    radeff_NO3 = 1E12 / 510072E9 * np.array([-166], dtype=dty)
elif (mod_NO3radeff == 'GEOS-CHEM'):
    radeff_NO3 = 1E12 / 510072E9 * np.array([-136], dtype=dty)
elif (mod_NO3radeff == 'GISS-MATRIX'):
    radeff_NO3 = 1E12 / 510072E9 * np.array([-240], dtype=dty)
elif (mod_NO3radeff == 'GMI'):
    radeff_NO3 = 1E12 / 510072E9 * np.array([-103], dtype=dty)
elif (mod_NO3radeff == 'HadGEM2'):
    radeff_NO3 = 1E12 / 510072E9 * np.array([-249], dtype=dty)
elif (mod_NO3radeff == 'IMPACT-Umich'):
    radeff_NO3 = 1E12 / 510072E9 * np.array([-155], dtype=dty)
elif (mod_NO3radeff == 'INCA'):
    radeff_NO3 = 1E12 / 510072E9 * np.array([-110], dtype=dty)
elif (mod_NO3radeff == 'NCAR-CAM-35'):
    radeff_NO3 = 1E12 / 510072E9 * np.array([-91], dtype=dty)
elif (mod_NO3radeff == 'OsloCTM2'):
    radeff_NO3 = 1E12 / 510072E9 * np.array([-173], dtype=dty)

# radiative efficiency of secondary organic aerosols {{W/m2}/Tg}
# from AeroCom2 [Myhre et al., 2013] (table 7)
if (mod_SOAradeff == 'mean-AeroCom2'):
    radeff_SOA = 1E12 / 510072E9 * np.array([-122], dtype=dty)
elif (mod_SOAradeff == 'CAM-51'):
    radeff_SOA = 1E12 / 510072E9 * np.array([-45], dtype=dty)
elif (mod_SOAradeff == 'GEOS-CHEM'):
    radeff_SOA = 1E12 / 510072E9 * np.array([-45], dtype=dty)
elif (mod_SOAradeff == 'IMPACT-Umich'):
    radeff_SOA = 1E12 / 510072E9 * np.array([-218], dtype=dty)
elif (mod_SOAradeff == 'MPIHAM'):
    radeff_SOA = 1E12 / 510072E9 * np.array([-139], dtype=dty)
elif (mod_SOAradeff == 'OsloCTM2'):
    radeff_SOA = 1E12 / 510072E9 * np.array([-161], dtype=dty)

# radiative efficiency of natural aerosols {{W/m2}/Tg}
# set to zero in this version
radeff_DUST = 1E12 / 510072E9 * np.array([0], dtype=dty)
radeff_SALT = 1E12 / 510072E9 * np.array([0], dtype=dty)

# ---------------
# 7.3.2. Indirect
# ---------------

# semi-direct effect (adjustements induced by the direct RF of BC)
# best-guess from IPCC AR5 [Boucher et al., 2013]
if (mod_BCadjust == 'Boucher2013'):
    k_BC_adjust = np.array([-0.1 / 0.6], dtype=dty)
# variations from [Lohmann et al., 2010] (figure 2; data provided by author)
elif (mod_BCadjust == 'CSIRO'):
    k_BC_adjust = np.array([-0.37], dtype=dty) * (-0.1 / 0.6) / -0.111
elif (mod_BCadjust == 'GISS'):
    k_BC_adjust = np.array([-0.225], dtype=dty) * (-0.1 / 0.6) / -0.111
elif (mod_BCadjust == 'HadGEM2'):
    k_BC_adjust = np.array([-0.13], dtype=dty) * (-0.1 / 0.6) / -0.111
elif (mod_BCadjust == 'ECHAM5'):
    k_BC_adjust = np.array([0.05], dtype=dty) * (-0.1 / 0.6) / -0.111
elif (mod_BCadjust == 'ECMWF'):
    k_BC_adjust = np.array([0.12], dtype=dty) * (-0.1 / 0.6) / -0.111

# solubility of aerosols for the aerosol-cloud interaction
# from [Hansen et al., 2005]
if (mod_CLOUDsolub == 'Hansen2005'):
    solub_SO4 = np.array([1.], dtype=dty)
    solub_POA = np.array([0.8], dtype=dty)
    solub_BC = np.array([0.7], dtype=dty)  # average of FF (0.6) and BB (0.8)
    solub_NO3 = np.array([1.], dtype=dty)
    solub_SOA = np.array([0.8], dtype=dty)  # assumed
    solub_DUST = np.array([0.], dtype=dty)
    solub_SALT = np.array([1.], dtype=dty)
# from [Lamarque et al., 2011] (RCP database)
elif (mod_CLOUDsolub == 'Lamarque2011'):
    solub_SO4 = np.array([1.], dtype=dty)
    solub_POA = np.array([0.86], dtype=dty)
    solub_BC = np.array([0.80], dtype=dty)
    solub_NO3 = np.array([1.], dtype=dty)
    solub_SOA = np.array([1.], dtype=dty)
    solub_DUST = np.array([0.12], dtype=dty)
    solub_SALT = np.array([0.05], dtype=dty)

# ERF over 1850-2000 for the aerosol-cloud interaction
# from ACCMIP [Shindell et al., 2013] (table 7)
# rescaled to the best guess of IPCC AR5 [Boucher et al., 2013]
if (mod_CLOUDerf == 'mean-ACCMIP'):
    RF_ref1 = np.array([-0.84], dtype=dty) * -0.45 / -0.84
elif (mod_CLOUDerf == 'CSIRO-Mk360'):
    RF_ref1 = np.array([-0.99], dtype=dty) * -0.45 / -0.84
elif (mod_CLOUDerf == 'GFDL-AM3'):
    RF_ref1 = np.array([-0.82], dtype=dty) * -0.45 / -0.84
elif (mod_CLOUDerf == 'GISS-E2-R'):
    RF_ref1 = np.array([-0.61], dtype=dty) * -0.45 / -0.84
elif (mod_CLOUDerf == 'HadGEM2'):
    RF_ref1 = np.array([-0.89], dtype=dty) * -0.45 / -0.84
elif (mod_CLOUDerf == 'LMDzORINCA'):
    RF_ref1 = np.array([-0.21], dtype=dty) * -0.45 / -0.84
elif (mod_CLOUDerf == 'MIROC-CHEM'):
    RF_ref1 = np.array([-1.12], dtype=dty) * -0.45 / -0.84
elif (mod_CLOUDerf == 'NCAR-CAM-51'):
    RF_ref1 = np.array([-1.22], dtype=dty) * -0.45 / -0.84

# soluble aerosol load for the aerosol-cloud interaction
# load pre-processed ACCMIP data
TMP = np.array([line for line in csv.reader(open(
    'data/AeroCloud_ACCMIP/#DATA.AeroCloud_' + mod_CLOUDerf + '.(2yr)_(7aer).LOAD.csv',
    'r'))][1:], dtype=dty)[:, 1:]
lgd = [line for line in csv.reader(open(
    'data/AeroCloud_ACCMIP/#DATA.AeroCloud_' + mod_CLOUDerf + '.(2yr)_(7aer).LOAD.csv',
    'r'))][0][1:]
AER_ref0 = 0
AER_ref1 = 0
for n in range(len(lgd)):
    if not np.isnan(TMP[0, n]):
        exec('AER_ref0 += TMP[0,n] * solub_' + lgd[n])
        exec('AER_ref1 += TMP[1,n] * solub_' + lgd[n])

# calculate parameters for aerosol-cloud interaction
# intensity of indirect effect {W/m2}
# based on a logarithmic formulation [e.g. Gultepe and Isaac, 1999]
Phi_0 = RF_ref1 / np.log(AER_ref1 / AER_ref0)

# preindustrial load of soluble aerosols {Tg}
# reduced by a factor from [Carslaw et al., 2013] and two arbitrary variations
if (mod_CLOUDpreind == 'median'):
    AERh_0 = AER_ref0 * np.exp(1 * (1.42 - 1.30) / Phi_0)
elif (mod_CLOUDpreind == 'high'):
    AERh_0 = AER_ref0 * np.exp(0 * (1.42 - 1.30) / Phi_0)
elif (mod_CLOUDpreind == 'low'):
    AERh_0 = AER_ref0 * np.exp(2 * (1.42 - 1.30) / Phi_0)

# ----------------
# 7.3.3. Volcanoes
# ----------------

# warming efficacy of volcano forcing {.}
# based on [Gregory et al., 2016]
warmeff_volc = np.array([0.6], dtype=dty)

# ===========
# 7.4. ALBEDO
# ===========

# -------------------
# 7.4.1. Black Carbon
# -------------------

# read region distribution
TMP = np.array([line for line in csv.reader(
    open('data/RegDiv_Reddy2007/#DATA.RegDiv_Reddy2007.114reg1_(9reg0).AREA.csv', 'r'))][
               1:], dtype=dty)
p_reg9 = np.zeros([nb_regionI, 9 + 1], dtype=dty)
for i in range(1, 114 + 1):
    p_reg9[regionI_index[i], :] += TMP[i - 1, :]
p_reg9 /= np.sum(p_reg9, 1)[:, np.newaxis]
p_reg9[np.isnan(p_reg9) | np.isinf(p_reg9)] = 0

# regional effect coefficient {.}
# from [Reddy and Boucher, 2007] (table 1)
if (mod_ALBBCreg == 'Reddy2007'):
    w_reg_BCsnow = np.array(
        [1, 1 / 6., 11 / 11., 1 / 10., 63 / 12., 2 / 3., 2 / 13., 17 / 43., 1 / 1.,
         2 / 1.], dtype=dty)

# global normalized radiative effect {{W/m2}/{Tg/yr}}
# from ACCMIP [Lee et al., 2013] (tab. 3 & fig. 15)
# rescaled to the best guess of IPCC AR5 [Boucher et al., 2013] (using table 7.1a)
if (mod_ALBBCrf == 'mean-ACCMIP'):
    radeff_BCsnow = np.array([0.0146 / (7.9 - 3.2)], dtype=dty) * (0.04 / 4.8) / (
                0.0146 / (7.9 - 3.2))
elif (mod_ALBBCrf == 'CICERO-OsloCTM2'):
    radeff_BCsnow = np.array([0.0131 / (7.8 - 3.1)], dtype=dty) * (0.04 / 4.8) / (
                0.0146 / (7.9 - 3.2))
elif (mod_ALBBCrf == 'GFDL-AM3'):
    radeff_BCsnow = np.array([0.0130 / (7.8 - 3.1)], dtype=dty) * (0.04 / 4.8) / (
                0.0146 / (7.9 - 3.2))
elif (mod_ALBBCrf == 'GISS-E2-R'):
    radeff_BCsnow = np.array([0.0142 / (8.8 - 4.0)], dtype=dty) * (0.04 / 4.8) / (
                0.0146 / (7.9 - 3.2))
elif (mod_ALBBCrf == 'GISS-E2-R-TOMAS'):
    radeff_BCsnow = np.array([0.0175 / (7.8 - 3.1)], dtype=dty) * (0.04 / 4.8) / (
                0.0146 / (7.9 - 3.2))
elif (mod_ALBBCrf == 'HadGEM2'):
    radeff_BCsnow = np.array([0.0133 / (7.8 - 3.1)], dtype=dty) * (0.04 / 4.8) / (
                0.0146 / (7.9 - 3.2))
elif (mod_ALBBCrf == 'MIROC-CHEM'):
    radeff_BCsnow = np.array([0.0173 / (7.7 - 3.0)], dtype=dty) * (0.04 / 4.8) / (
                0.0146 / (7.9 - 3.2))
elif (mod_ALBBCrf == 'NCAR-CAM-35'):
    radeff_BCsnow = np.array([0.0143 / (7.8 - 3.1)], dtype=dty) * (0.04 / 4.8) / (
                0.0146 / (7.9 - 3.2))
elif (mod_ALBBCrf == 'NCAR-CAM-51'):
    radeff_BCsnow = np.array([0.0141 / (7.8 - 3.1)], dtype=dty) * (0.04 / 4.8) / (
                0.0146 / (7.9 - 3.2))

# warming efficacy of black carbon on snow {.}
# from IPCC AR5 [Boucher et al., 2013] (sect. 7.5.2.3)
if (mod_ALBBCwarm == 'median'):
    warmeff_BCsnow = np.array([3.], dtype=dty)
elif (mod_ALBBCwarm == 'low'):
    warmeff_BCsnow = np.array([2.], dtype=dty)
elif (mod_ALBBCwarm == 'high'):
    warmeff_BCsnow = np.array([4.], dtype=dty)

# -----------------
# 7.4.2. Land-Cover
# -----------------

# basic biomes of aggregation
bio = ['des', 'for', 'shr', 'gra', 'cro', 'pas', 'urb']

# load pre-processed albedo climatology
TMP = np.array([line for line in csv.reader(open(
    'data/Albedo_' + mod_ALBLCalb + '/#DATA.Albedo_' + mod_ALBLCalb + '_' + mod_ALBLCflux + '_' + mod_ALBLCcover + '.114reg1_7bio.alb.csv',
    'r'))], dtype=dty)
TMP2 = np.array([line for line in csv.reader(open(
    'data/Albedo_' + mod_ALBLCalb + '/#DATA.Albedo_' + mod_ALBLCalb + '_' + mod_ALBLCflux + '_' + mod_ALBLCcover + '.114reg1_7bio.RSDS.csv',
    'r'))], dtype=dty)
alpha_alb = np.zeros([nb_regionI, nb_biome])
RSDS_alb = np.zeros([nb_regionI, nb_biome])
for i in range(1, 114 + 1):
    for b in range(len(bio)):
        alpha_alb[regionI_index[i], biome_index[bio[b]]] += TMP[i - 1, b] * TMP2[i - 1, b]
        RSDS_alb[regionI_index[i], biome_index[bio[b]]] += TMP2[i - 1, b]
alpha_alb /= RSDS_alb
alpha_alb[np.isnan(alpha_alb) | np.isinf(alpha_alb)] = 0

# set pasture albedo
if (np.sum(alpha_alb[:, biome_index["pas"]]) == 0):
    for i in range(1, 114 + 1):
        alpha_alb[regionI_index[i], biome_index["pas"]] += 0.6 * TMP[
            i - 1, bio.index("gra")] * TMP2[i - 1, bio.index("gra")] + 0.4 * TMP[
                                                               i - 1, bio.index("des")] * \
                                                           TMP2[i - 1, bio.index("des")]
        RSDS_alb[regionI_index[i], biome_index["pas"]] += 0.6 * TMP2[
            i - 1, bio.index("gra")] + 0.4 * TMP2[i - 1, bio.index("des")]
    alpha_alb[:, biome_index["pas"]] /= RSDS_alb[:, biome_index["pas"]]
    alpha_alb[np.isnan(alpha_alb) | np.isinf(alpha_alb)] = 0

# load pre-processed radiation climatology
TMP = np.array([line for line in csv.reader(open(
    'data/RadFlux_' + mod_ALBLCflux + '/#DATA.RadFlux_' + mod_ALBLCflux + '.114reg1.rsds.csv',
    'r'))], dtype=dty)
TMP2 = np.array([line for line in csv.reader(open(
    'data/RadFlux_' + mod_ALBLCflux + '/#DATA.RadFlux_' + mod_ALBLCflux + '.114reg1.AREA.csv',
    'r'))], dtype=dty)
rsds_alb = np.zeros([nb_regionI])
AREA_alb = np.zeros([nb_regionI])
for i in range(1, 114 + 1):
    rsds_alb[regionI_index[i]] += TMP[i - 1, 0] * TMP2[i - 1, 0]
    AREA_alb[regionI_index[i]] += TMP2[i - 1, 0]
rsds_alb /= AREA_alb
rsds_alb[np.isnan(rsds_alb) | np.isinf(rsds_alb)] = 0

# upward transmittance from [Lenton and Vaughan, 2009]
p_trans = np.array([-0.854], dtype=dty)

# final albedo parameters {{W/m2}/Mha}
alpha_LCC = p_trans * alpha_alb * rsds_alb[:, np.newaxis] / (510072E9 / 1E10)

# warming efficacy of land-cover change albedo effect {.}
# from [Bright et al., 2015] (tab. 7)
if (mod_ALBLCwarm == 'Hansen2005'):
    warmeff_LCC = np.array([1.02], dtype=dty)
elif (mod_ALBLCwarm == 'Davin2007'):
    warmeff_LCC = np.array([0.5], dtype=dty)
elif (mod_ALBLCwarm == 'Davin2010'):
    warmeff_LCC = np.array([0.78], dtype=dty)
elif (mod_ALBLCwarm == 'Jones2013'):
    warmeff_LCC = np.array([0.79], dtype=dty)

##################################################
#   8. CLIMATE
##################################################

# ==========
# 8.1. GLOBE
# ==========

# ----------------------
# 8.1.A. Reconstructions
# ----------------------

# historical global temperatures from NOAA/NCDC {degC}
# from [Smith et al., 2008] updated from website
for var in ['gst', 'lst', 'sst']:
    exec(var + '_ncdc = np.ones([314+1], dtype=dty) * np.nan')
    TMP = np.array([line for line in csv.reader(
        open('data/HistClim_NOAA-NCDC/#DATA.HistClim_NOAA-NCDC.1880-2014.' + var + '.csv',
             'r'))], dtype=dty)
    exec(var + '_ncdc[180:] = TMP[:,0]')

# historical global temperature from GISTEMP {degC}
# from [Hansen et al., 2010] updated from website
for var in ['gst']:
    exec(var + '_giss = np.ones([314+1], dtype=dty) * np.nan')
    TMP = np.array([line for line in csv.reader(
        open('data/HistClim_GISTEMP/#DATA.HistClim_GISTEMP.1880-2014.' + var + '.csv',
             'r'))], dtype=dty)
    exec(var + '_giss[180:] = TMP[:,0]')

# historical global temperature from HadCRUT4 {degC}
# from [Morice et al., 2012] updated from website
for var in ['gst']:
    exec(var + '_had = np.ones([314+1], dtype=dty) * np.nan')
    TMP = np.array([line for line in csv.reader(
        open('data/HistClim_HadCRUT4/#DATA.HistClim_HadCRUT4.1850-2014.' + var + '.csv',
             'r'))], dtype=dty)
    exec(var + '_had[150:] = TMP[:,0]')

# ------------
# 8.1.1. CMIP5
# ------------

# length of CMIP5 simulations per model
lng = {'mean-CMIP5': 140, 'ACCESS-10': 150, 'ACCESS-13': 151, 'BCC-CSM-11': 150,
       'BCC-CSM-11m': 150, 'CanESM2': 150, 'CCSM4': 151, 'CNRM-CM5': 150,
       'CNRM-CM5-2': 140, 'CSIRO-Mk360': 150, 'GFDL-CM3': 150, 'GFDL-ESM2G': 300,
       'GFDL-ESM2M': 300, 'GISS-E2-H': 151, 'GISS-E2-R': 151, 'HadGEM2-ES': 151,
       'IPSL-CM5A-LR': 260, 'IPSL-CM5A-MR': 140, 'IPSL-CM5B-LR': 160, 'MIROC5': 151,
       'MIROC-ESM': 150, 'MPI-ESM-LR': 150, 'MPI-ESM-MR': 150, 'MPI-ESM-P': 150,
       'MRI-CGCM3': 150, 'NorESM1-M': 150}

# load pre-processed CMIP5 results for specified model
# data related to temperature change
for sim in ['ctrl', 'quad']:
    TMP = np.array([line for line in csv.reader(open(
        'data/Climate_CMIP5/#DATA.Climate_' + mod_TEMPresp + '.' + str(
            lng[mod_TEMPresp]) + 'yr_(7var).' + sim + '_global.csv', 'r'))][1:],
                   dtype=dty)
    lgd = [line for line in csv.reader(open(
        'data/Climate_CMIP5/#DATA.Climate_' + mod_TEMPresp + '.' + str(
            lng[mod_TEMPresp]) + 'yr_(7var).' + sim + '_global.csv', 'r'))][0]
    exec('gst_' + sim + 'T = TMP[:,lgd.index("tas")]')
    exec(
        'erb_' + sim + 'T = TMP[:,lgd.index("rsdt")] - TMP[:,lgd.index("rsut")] - TMP[:,lgd.index("rlut")]')
# data related to precipitations change
for sim in ['ctrl', 'quad']:
    TMP = np.array([line for line in csv.reader(open(
        'data/Climate_CMIP5/#DATA.Climate_' + mod_PRECresp + '.' + str(
            lng[mod_PRECresp]) + 'yr_(7var).' + sim + '_global.csv', 'r'))][1:],
                   dtype=dty)
    lgd = [line for line in csv.reader(open(
        'data/Climate_CMIP5/#DATA.Climate_' + mod_PRECresp + '.' + str(
            lng[mod_PRECresp]) + 'yr_(7var).' + sim + '_global.csv', 'r'))][0]
    exec('gst_' + sim + 'P = TMP[:,lgd.index("tas")]')
    exec('gyp_' + sim + 'P = TMP[:,lgd.index("pr")]')

# definition of parameters
# equilibrium climate sensitivity {K}&{K/{W/m2}}
lambda_0 = np.array([0], dtype=dty)
# dynamical parameters for temperature change {.}&{yr}&{yr}
theta_0 = np.array([0], dtype=dty)
tau_gst = np.array([0], dtype=dty)
tau_gst0 = np.array([0], dtype=dty)
# precipitation sensitivities to temperature and radiative forcing {mm/K}&{mm/{W/m2}}
alpha_gyp = np.array([0], dtype=dty)
beta_gyp = np.array([0], dtype=dty)

# calculation of parameters
# equilibrium temperature change
# following the approach by [Gregory et al., 2004]
diff = gst_quadT - np.mean(gst_ctrlT, 0)


def err(var):
    clim = var[0] * (erb_quadT - np.mean(erb_ctrlT, 0)) + var[1]
    return np.sum((diff - clim) ** 2)


[TMP, T4x] = fmin(err, [-1., 5.], disp=False)
lambda_0[0] = T4x / f_RF_CO2(3 * CO2_0)

# impulse response function for temperature change
diff = gst_quadT - np.mean(gst_ctrlT, 0)


def err(var):
    clim = T4x * (1 - var[2] * np.exp(
        -np.arange(0.5, lng[mod_TEMPresp] + 0.5, 1) / var[0]) - (1 - var[2]) * np.exp(
        -np.arange(0.5, lng[mod_TEMPresp] + 0.5, 1) / var[1]))
    return np.sum((diff - clim) ** 2)


[tau_Tfast, tau_Tslow, p_Tfast] = fmin(err, [4., 400., 0.5], disp=False)
p_Tslow = 1 - p_Tfast
if tau_Tfast > tau_Tslow:
    TMP = tau_Tfast
    tau_Tfast = tau_Tslow
    tau_Tslow = TMP
    TMP = p_Tfast
    p_Tfast = p_Tslow
    p_Tslow = TMP

# parameters of corresponding two-box model
# formulation based on [Geoffroy et al., 2013]
tau_gst[0] = tau_Tfast * tau_Tslow / (tau_Tfast * p_Tslow + tau_Tslow * p_Tfast)
tau_gst0[0] = tau_Tfast * p_Tfast + tau_Tslow * p_Tslow - tau_gst[0]
theta_0[0] = tau_gst0[0] / (tau_Tfast * p_Tslow + tau_Tslow * p_Tfast)

# sensitivities of precipitation change
diff = gyp_quadP - np.mean(gyp_ctrlP, 0)


def err(var):
    clim = np.abs(var[0]) * (gst_quadP - np.mean(gst_ctrlP, 0)) - np.abs(var[1])
    return np.sum((diff - clim) ** 2)


[alpha_gyp[0], beta_gyp[0]] = fmin(err, [10., -50.], disp=False)
alpha_gyp[0] = np.abs(alpha_gyp[0])
beta_gyp[0] = -np.abs(beta_gyp[0]) / f_RF_CO2(3 * CO2_0)

# ---------------------
# 8.1.2. Precipitations
# ---------------------

# global precipitation response differentiated by forcing agent
# formulation based on [Allan et al., 2014]

# factors from [Andrews et al., 2010] (tab. 3)
# assumptions following [Allan et al., 2014]
if (mod_PRECradfact == 'Andrews2010'):
    p_atm_CO2 = np.array([0.8], dtype=dty)
    p_atm_noCO2 = np.array([0.5], dtype=dty)
    p_atm_O3t = np.array([-0.3], dtype=dty)
    p_atm_strat = np.array([0.0], dtype=dty)  # assumed
    p_atm_scatter = np.array([0.0], dtype=dty)
    p_atm_absorb = np.array([2.5], dtype=dty)
    p_atm_cloud = np.array([0.0], dtype=dty)  # assumed
    p_atm_alb = np.array([0.0], dtype=dty)
    p_atm_solar = np.array([0.2], dtype=dty)

# factors from [Kvalevag et al., 2013] (tab. 2; highest perturbation)
# assumptions following [Allan et al., 2014] and nil if not given
elif (mod_PRECradfact == 'Kvalevag2013'):
    p_atm_CO2 = np.array([2.0 / 3.6], dtype=dty)
    p_atm_noCO2 = np.array([0.5 / 1.8], dtype=dty)
    p_atm_O3t = np.array([0.0], dtype=dty)  # not given
    p_atm_strat = np.array([0.0], dtype=dty)  # assumed
    p_atm_scatter = np.array([0.6 / -1.4], dtype=dty)
    p_atm_absorb = np.array([3.1 / 0.5], dtype=dty)
    p_atm_cloud = np.array([0.0], dtype=dty)  # assumed
    p_atm_alb = np.array([0.0], dtype=dty)  # not given
    p_atm_solar = np.array([0.4 / 3.5], dtype=dty)

# normalization of radiative factor of precipitations
beta_gyp /= p_atm_CO2

# ==========
# 8.2. OCEAN
# ==========

# ----------------------
# 8.2.A. Reconstructions
# ----------------------

# historical sea climate from HadISST1 {degC}&{Mha}
# from [Rayner et al., 2003] updated from website
for var in ['sst', 'sic']:
    exec(var + '_had = np.ones([314+1], dtype=dty) * np.nan')
    TMP = np.array([line for line in csv.reader(open(
        'data/HistOcean_HadISST1/#DATA.HistOcean_HadISST1.1870-2014_18x10lat.' + var + '.csv',
        'r'))], dtype=dty)
    TMP2 = np.array([line for line in csv.reader(open(
        'data/HistOcean_HadISST1/#DATA.HistOcean_HadISST1.1870-2014_18x10lat.SURF.csv',
        'r'))], dtype=dty)
    exec(var + '_had[170:] = np.sum(TMP[:,:]*TMP2[:,:],1)')
    exec(var + '_had[170:] /= np.sum(TMP2[:,:],1)')

# historical sea climate from ERSST4 {degC}
# from [Huang et al., 2015] updated from website
for var in ['sst']:
    exec(var + '_ncdc2 = np.ones([314+1], dtype=dty) * np.nan')
    TMP = np.array([line for line in csv.reader(open(
        'data/HistOcean_ERSST4/#DATA.HistOcean_ERSST4.1854-2014_18x10lat.' + var + '.csv',
        'r'))], dtype=dty)
    TMP2 = np.array([line for line in csv.reader(
        open('data/HistOcean_ERSST4/#DATA.HistOcean_ERSST4.1854-2014_18x10lat.SURF.csv',
             'r'))], dtype=dty)
    exec(var + '_ncdc2[154:] = np.sum(TMP[:,:]*TMP2[:,:],1)')
    exec(var + '_ncdc2[154:] /= np.sum(TMP2[:,:],1)')

# historical ocean heat content
# from [Levitus et al., 2012] updated from website
# reference period is whole period
for var in ['D_OHC']:
    exec(var + '_nodc = np.ones([311+1], dtype=dty) * np.nan')
    TMP = np.array([line for line in csv.reader(
        open('data/HistOcean_NOAA-NODC/#DATA.HistOcean_NOAA-NODC.1955-2011.D_OHC.csv',
             'r'))], dtype=dty)
    exec(var + '_nodc[255:] = TMP[:,0]')

# ------------
# 8.2.1. CMIP5
# ------------

# length of CMIP5 simulations per model or simulation
lng = {'mean-CMIP5': 140, 'ACCESS-10': 150, 'ACCESS-13': 151, 'BCC-CSM-11': 150,
       'BCC-CSM-11m': 150, 'CanESM2': 150, 'CCSM4': 151, 'CNRM-CM5': 150,
       'CNRM-CM5-2': 140, 'CSIRO-Mk360': 150, 'GFDL-CM3': 150, 'GFDL-ESM2G': 300,
       'GFDL-ESM2M': 300, 'GISS-E2-H': 151, 'GISS-E2-R': 151, 'HadGEM2-ES': 151,
       'IPSL-CM5A-LR': 260, 'IPSL-CM5A-MR': 140, 'IPSL-CM5B-LR': 160, 'MIROC5': 151,
       'MIROC-ESM': 150, 'MPI-ESM-LR': 150, 'MPI-ESM-MR': 150, 'MPI-ESM-P': 150,
       'MRI-CGCM3': 150, 'NorESM1-M': 150}
prd = {'ctrl': '251yr', 'hist': '1850-2005', 'rcp26': '2006-2100', 'rcp45': '2006-2100',
       'rcp60': '2006-2100', 'rcp85': '2006-2100'}

# load pre-processed CMIP5 results for specified model
# temperature patterns based on 'abrupt4xCO2'
if (mod_TEMPpattern == '4xCO2'):
    for sim in ['ctrl', 'quad']:
        # global
        TMP = np.array([line for line in csv.reader(open(
            'data/Climate_CMIP5/#DATA.Climate_' + mod_TEMPresp + '.' + str(
                lng[mod_TEMPresp]) + 'yr_(7var).' + sim + '_global.csv', 'r'))][1:],
                       dtype=dty)
        lgd = [line for line in csv.reader(open(
            'data/Climate_CMIP5/#DATA.Climate_' + mod_TEMPresp + '.' + str(
                lng[mod_TEMPresp]) + 'yr_(7var).' + sim + '_global.csv', 'r'))][0]
        exec('gst_' + sim + 'TR = TMP[:,lgd.index("tas")]')
        # local
        for var in ['tos']:
            TMP = np.array([line for line in csv.reader(open(
                'data/Climate_CMIP5/#DATA.Climate_' + mod_TEMPresp + '.' + str(
                    lng[mod_TEMPresp]) + 'yr_18x10lat.' + sim + '_' + var + '.csv',
                'r'))], dtype=dty)
            TMP2 = np.array([line for line in csv.reader(open(
                'data/Climate_CMIP5/#DATA.Climate_' + mod_TEMPresp + '.' + str(
                    lng[mod_TEMPresp]) + 'yr_18x10lat.SURF.csv', 'r'))], dtype=dty)
            exec(var + '_' + sim + 'TR = np.sum(TMP*TMP2,1)/np.sum(TMP2,1)')

# temperature patterns based on 'historical' and 'rcp'
elif (mod_TEMPpattern == 'hist&RCPs'):
    for sim in ['ctrl', 'hist', 'rcp26', 'rcp45', 'rcp60', 'rcp85']:
        # global
        if os.path.isfile('data/ClimReg_CMIP5/#DATA.ClimReg_' + mod_TEMPresp + '.' + prd[
            sim] + '_(3var).' + sim + '_global.csv'):
            TMP = np.array([line for line in csv.reader(open(
                'data/ClimReg_CMIP5/#DATA.ClimReg_' + mod_TEMPresp + '.' + prd[
                    sim] + '_(3var).' + sim + '_global.csv', 'r'))][1:], dtype=dty)
            lgd = [line for line in csv.reader(open(
                'data/ClimReg_CMIP5/#DATA.ClimReg_' + mod_TEMPresp + '.' + prd[
                    sim] + '_(3var).' + sim + '_global.csv', 'r'))][0]
            exec('gst_' + sim + 'TR = TMP[:,lgd.index("tas")]')
        # local
        for var in ['tos']:
            if os.path.isfile(
                    'data/ClimReg_CMIP5/#DATA.ClimReg_' + mod_TEMPresp + '.' + prd[
                        sim] + '_18x10lat.' + sim + '_' + var + '.csv'):
                TMP = np.array([line for line in csv.reader(open(
                    'data/ClimReg_CMIP5/#DATA.ClimReg_' + mod_TEMPresp + '.' + prd[
                        sim] + '_18x10lat.' + sim + '_' + var + '.csv', 'r'))], dtype=dty)
                TMP2 = np.array([line for line in csv.reader(open(
                    'data/ClimReg_CMIP5/#DATA.ClimReg_' + mod_TEMPresp + '.251yr_18x10lat.SURF.csv',
                    'r'))], dtype=dty)
                exec(
                    var + '_' + sim + 'TR = np.sum(TMP*TMP2[:len(TMP)],1)/np.sum(TMP2[:len(TMP)],1)')

# aggregate all experiments
for VAR in ['gst', 'tos']:
    if (mod_TEMPpattern == '4xCO2'):
        exec(VAR + '_allTR = np.array(list(' + VAR + '_quadTR), dtype=dty)')
    elif (mod_TEMPpattern == 'hist&RCPs'):
        exec(VAR + '_allTR = []')
        for sim in ['rcp26', 'rcp45', 'rcp60', 'rcp85']:
            if os.path.isfile(
                    'data/ClimReg_CMIP5/#DATA.ClimReg_' + mod_TEMPresp + '.' + prd[
                        sim] + '_(3var).' + sim + '_global.csv'):
                exec(
                    VAR + '_allTR += list(' + VAR + '_histTR)+list(' + VAR + '_' + sim + 'TR)')
        exec('test = ' + VAR + '_allTR')
        if (test == []):
            exec(VAR + '_allTR += list(' + VAR + '_histTR)')
        exec(VAR + '_allTR = np.array(' + VAR + '_allTR, dtype=dty)')

# definition of parameter
# scaling of sea surface temperature over gst {.}
w_reg_sst = np.array([0], dtype=dty)

# fit of parameter
diff = tos_allTR - np.mean(tos_ctrlTR, 0)


def err(var):
    clim = var[0] * (gst_allTR - np.mean(gst_ctrlTR, 0))
    return np.sum((diff - clim) ** 2)


[w_reg_sst[0]] = fmin(err, [1], disp=False)

# -------------------
# 8.2.2. Heat content
# -------------------

# conversion factor from radiative forcing {W/m2} to heat {ZJ}
alpha_OHC = np.array([510072E9 * 3600 * 24 * 365.25 / 1E21], dtype=dty)

# fraction of energy used to heat the land, the atmosphere, and to melt ice {.}
# from [Otto et al., 2013] (SI)
p_OHC = np.array([1 - 0.03 - 0.01 - 0.02], dtype=dty)

# --------------------
# 8.2.3. Acidification
# --------------------

# fonction relating surface pH and atmospheric CO2
# from [Tans, 2009] based on [Tans, 2008]
if (mod_ACIDsurf == 'Tans2009'):

    def f_pH(D_CO2):
        D_pH = -0.85 * np.log(1 + D_CO2 / CO2_0)
        return np.array(D_pH, dtype=dty)

# from [Bernie et al., 2010]
elif (mod_ACIDsurf == 'Bernie2010'):

    def f_pH(D_CO2):
        D_pH = D_pH = -0.00173 * D_CO2 + 1.3264E-6 * (
                    2 * CO2_0 * D_CO2 + D_CO2 ** 2) - 4.4943E-10 * (
                                  3 * D_CO2 * CO2_0 ** 2 + 3 * CO2_0 * D_CO2 ** 2 + D_CO2 ** 3)
        return np.array(D_pH, dtype=dty)

# ----------------
# 8.2.4. Sea-level
# ----------------

# TODO

# =========
# 8.3. LAND
# =========

# ----------------------
# 8.3.A. Reconstructions
# ----------------------

# historical local climate from CRU {degC}&{mm}
# from [Harris et al., 2014]
for var in ['lst', 'lyp']:
    exec(var + '_cru = np.zeros([314+1,nb_regionI], dtype=dty)')
    TMP = np.array([line for line in csv.reader(
        open('data/HistLand_CRU/#DATA.HistLand_CRU.1901-2014_114reg1.' + var + '.csv',
             'r'))], dtype=dty)
    TMP2 = np.array([line for line in csv.reader(
        open('data/HistLand_CRU/#DATA.HistLand_CRU.1901-2014_114reg1.AREA.csv', 'r'))],
                    dtype=dty)
    for i in range(1, 114 + 1):
        exec(var + '_cru[201:,regionI_index[i]] += TMP[:,i-1]*TMP2[:,i-1]')
    TMP = np.zeros([314 + 1, nb_regionI], dtype=dty)
    for i in range(1, 114 + 1):
        TMP[201:, regionI_index[i]] += TMP2[:, i - 1]
    exec(var + '_cru /= TMP[:,:]')

# historical local climate from GHCN+CAMS {degC}
# from [Fan and van den Dool, 2008] updated from website
for var in ['lst']:
    exec(var + '_ghcn = np.zeros([314+1,nb_regionI], dtype=dty)')
    TMP = np.array([line for line in csv.reader(open(
        'data/HistLand_GHCN-CAMS/#DATA.HistLand_GHCN-CAMS.1948-2014_114reg1.' + var + '.csv',
        'r'))], dtype=dty)
    TMP2 = np.array([line for line in csv.reader(open(
        'data/HistLand_GHCN-CAMS/#DATA.HistLand_GHCN-CAMS.1948-2014_114reg1.AREA.csv',
        'r'))], dtype=dty)
    for i in range(1, 114 + 1):
        exec(var + '_ghcn[248:,regionI_index[i]] += TMP[:,i-1]*TMP2[:,i-1]')
    TMP = np.zeros([314 + 1, nb_regionI], dtype=dty)
    for i in range(1, 114 + 1):
        TMP[248:, regionI_index[i]] += TMP2[:, i - 1]
    exec(var + '_ghcn /= TMP[:,:]')

# ------------
# 8.3.1. CMIP5
# ------------

# length of CMIP5 simulations per model or simulation
lng = {'mean-CMIP5': 140, 'ACCESS-10': 150, 'ACCESS-13': 151, 'BCC-CSM-11': 150,
       'BCC-CSM-11m': 150, 'CanESM2': 150, 'CCSM4': 151, 'CNRM-CM5': 150,
       'CNRM-CM5-2': 140, 'CSIRO-Mk360': 150, 'GFDL-CM3': 150, 'GFDL-ESM2G': 300,
       'GFDL-ESM2M': 300, 'GISS-E2-H': 151, 'GISS-E2-R': 151, 'HadGEM2-ES': 151,
       'IPSL-CM5A-LR': 260, 'IPSL-CM5A-MR': 140, 'IPSL-CM5B-LR': 160, 'MIROC5': 151,
       'MIROC-ESM': 150, 'MPI-ESM-LR': 150, 'MPI-ESM-MR': 150, 'MPI-ESM-P': 150,
       'MRI-CGCM3': 150, 'NorESM1-M': 150}
prd = {'ctrl': '251yr', 'hist': '1850-2005', 'rcp26': '2006-2100', 'rcp45': '2006-2100',
       'rcp60': '2006-2100', 'rcp85': '2006-2100'}

# load pre-processed CMIP5 results for specified model
# temperature patterns based on 'abrupt4xCO2'
if (mod_TEMPpattern == '4xCO2'):
    for sim in ['ctrl', 'quad']:
        # global
        TMP = np.array([line for line in csv.reader(open(
            'data/Climate_CMIP5/#DATA.Climate_' + mod_TEMPresp + '.' + str(
                lng[mod_TEMPresp]) + 'yr_(7var).' + sim + '_global.csv', 'r'))][1:],
                       dtype=dty)
        lgd = [line for line in csv.reader(open(
            'data/Climate_CMIP5/#DATA.Climate_' + mod_TEMPresp + '.' + str(
                lng[mod_TEMPresp]) + 'yr_(7var).' + sim + '_global.csv', 'r'))][0]
        exec('gst_' + sim + 'TR = TMP[:,lgd.index("tas")]')
        # local
        for var in ['tas']:
            exec(
                var + '_' + sim + 'TR = np.zeros([lng[mod_TEMPresp],nb_regionI], dtype=dty)')
            TMP = np.array([line for line in csv.reader(open(
                'data/Climate_CMIP5/#DATA.Climate_' + mod_TEMPresp + '.' + str(
                    lng[mod_TEMPresp]) + 'yr_114reg1.' + sim + '_' + var + '.csv', 'r'))],
                           dtype=dty)
            TMP2 = np.array([line for line in csv.reader(open(
                'data/Climate_CMIP5/#DATA.Climate_' + mod_TEMPresp + '.' + str(
                    lng[mod_TEMPresp]) + 'yr_114reg1.AREA.csv', 'r'))], dtype=dty)
            exec(var + '_' + sim + 'TR = np.zeros([len(TMP),nb_regionI], dtype=dty)')
            exec('AREA_' + sim + 'TR = np.zeros([len(TMP),nb_regionI], dtype=dty)')
            for i in range(1, 114 + 1):
                exec(
                    var + '_' + sim + 'TR[:,regionI_index[i]] += TMP[:,i-1]*TMP2[:len(TMP),i-1]')
                exec('AREA_' + sim + 'TR[:,regionI_index[i]] += TMP2[:len(TMP),i-1]')
            exec(var + '_' + sim + 'TR /= AREA_' + sim + 'TR')
            exec(
                var + '_' + sim + 'TR[np.isnan(' + var + '_' + sim + 'TR)|np.isinf(' + var + '_' + sim + 'TR)] = 0')

# temperature patterns based on 'historical' and 'rcp'
elif (mod_TEMPpattern == 'hist&RCPs'):
    for sim in ['ctrl', 'hist', 'rcp26', 'rcp45', 'rcp60', 'rcp85']:
        # global
        if os.path.isfile('data/ClimReg_CMIP5/#DATA.ClimReg_' + mod_TEMPresp + '.' + prd[
            sim] + '_(3var).' + sim + '_global.csv'):
            TMP = np.array([line for line in csv.reader(open(
                'data/ClimReg_CMIP5/#DATA.ClimReg_' + mod_TEMPresp + '.' + prd[
                    sim] + '_(3var).' + sim + '_global.csv', 'r'))][1:], dtype=dty)
            lgd = [line for line in csv.reader(open(
                'data/ClimReg_CMIP5/#DATA.ClimReg_' + mod_TEMPresp + '.' + prd[
                    sim] + '_(3var).' + sim + '_global.csv', 'r'))][0]
            exec('gst_' + sim + 'TR = TMP[:,lgd.index("tas")]')
        # local
        for var in ['tas']:
            if os.path.isfile(
                    'data/ClimReg_CMIP5/#DATA.ClimReg_' + mod_TEMPresp + '.' + prd[
                        sim] + '_114reg1.' + sim + '_' + var + '.csv'):
                TMP = np.array([line for line in csv.reader(open(
                    'data/ClimReg_CMIP5/#DATA.ClimReg_' + mod_TEMPresp + '.' + prd[
                        sim] + '_114reg1.' + sim + '_' + var + '.csv', 'r'))], dtype=dty)
                TMP2 = np.array([line for line in csv.reader(open(
                    'data/ClimReg_CMIP5/#DATA.ClimReg_' + mod_TEMPresp + '.251yr_114reg1.AREA.csv',
                    'r'))], dtype=dty)
                exec(var + '_' + sim + 'TR = np.zeros([len(TMP),nb_regionI], dtype=dty)')
                exec('AREA_' + sim + 'TR = np.zeros([len(TMP),nb_regionI], dtype=dty)')
                for i in range(1, 114 + 1):
                    exec(
                        var + '_' + sim + 'TR[:,regionI_index[i]] += TMP[:,i-1]*TMP2[:len(TMP),i-1]')
                    exec('AREA_' + sim + 'TR[:,regionI_index[i]] += TMP2[:len(TMP),i-1]')
                exec(var + '_' + sim + 'TR /= AREA_' + sim + 'TR')
                exec(
                    var + '_' + sim + 'TR[np.isnan(' + var + '_' + sim + 'TR)|np.isinf(' + var + '_' + sim + 'TR)] = 0')

# aggregate all experiments
for VAR in ['gst', 'tas']:
    if (mod_TEMPpattern == '4xCO2'):
        exec(VAR + '_allTR = np.array(list(' + VAR + '_quadTR), dtype=dty)')
    elif (mod_TEMPpattern == 'hist&RCPs'):
        exec(VAR + '_allTR = []')
        nrcp = 0
        for sim in ['rcp26', 'rcp45', 'rcp60', 'rcp85']:
            if os.path.isfile(
                    'data/ClimReg_CMIP5/#DATA.ClimReg_' + mod_TEMPresp + '.' + prd[
                        sim] + '_(3var).' + sim + '_global.csv'):
                nrcp += 1
                exec(
                    VAR + '_allTR += list(' + VAR + '_histTR)+list(' + VAR + '_' + sim + 'TR)')
        exec('test = ' + VAR + '_allTR')
        if (test == []):
            exec(VAR + '_allTR += list(' + VAR + '_histTR)')
        exec(VAR + '_allTR = np.array(' + VAR + '_allTR, dtype=dty)')
# decadal means
for VAR in ['gst', 'tas']:
    if (mod_TEMPpattern == '4xCO2'):
        exec(
            VAR + '_decTR= np.zeros([lng[mod_TEMPresp]-10]+list(np.shape(' + VAR + '_allTR)[1:]), dtype=dty)')
        for t in range(lng[mod_TEMPresp] - 10):
            exec(VAR + '_decTR[t,...] = np.mean(' + VAR + '_allTR[t:t+10,...],0)')
    elif (mod_TEMPpattern == 'hist&RCPs'):
        exec(
            VAR + '_decTR= np.zeros([(156-10)+nrcp*(95-10)]+list(np.shape(' + VAR + '_allTR)[1:]), dtype=dty)')
        for t in range(156 - 10):
            exec(VAR + '_decTR[t,...] = np.mean(' + VAR + '_allTR[t:t+10,...],0)')
        for n in range(nrcp):
            for t in range(95 - 10):
                exec(
                    VAR + '_decTR[(156-10)+n*(95-10)+t,...] = np.mean(' + VAR + '_allTR[(156-10)+n*(95-10)+t:(156-10)+n*(95-10)+t+10,...],0)')

            # definition of parameter
# scaling of local surface temperature over gst {.}
w_reg_lst = np.zeros([nb_regionI], dtype=dty)

# fit of parameter
for i in range(nb_regionI):
    diff = tas_decTR[:, i] - np.mean(tas_ctrlTR[:, i], 0)


    def err(var):
        clim = var[0] * (gst_decTR - np.mean(gst_ctrlTR, 0))
        return np.sum((diff - clim) ** 2)


    [w_reg_lst[i]] = fmin(err, [1], disp=False)

# load pre-processed CMIP5 results for specified model
# precipitation patterns based on 'abrupt4xCO2'
if (mod_PRECpattern == '4xCO2'):
    for sim in ['ctrl', 'quad']:
        # global
        TMP = np.array([line for line in csv.reader(open(
            'data/Climate_CMIP5/#DATA.Climate_' + mod_PRECresp + '.' + str(
                lng[mod_PRECresp]) + 'yr_(7var).' + sim + '_global.csv', 'r'))][1:],
                       dtype=dty)
        lgd = [line for line in csv.reader(open(
            'data/Climate_CMIP5/#DATA.Climate_' + mod_PRECresp + '.' + str(
                lng[mod_PRECresp]) + 'yr_(7var).' + sim + '_global.csv', 'r'))][0]
        exec('gyp_' + sim + 'PR = TMP[:,lgd.index("pr")]')
        # local
        for var in ['pr']:
            TMP = np.array([line for line in csv.reader(open(
                'data/Climate_CMIP5/#DATA.Climate_' + mod_PRECresp + '.' + str(
                    lng[mod_PRECresp]) + 'yr_114reg1.' + sim + '_' + var + '.csv', 'r'))],
                           dtype=dty)
            TMP2 = np.array([line for line in csv.reader(open(
                'data/Climate_CMIP5/#DATA.Climate_' + mod_PRECresp + '.' + str(
                    lng[mod_PRECresp]) + 'yr_114reg1.AREA.csv', 'r'))], dtype=dty)
            exec(var + '_' + sim + 'PR = np.zeros([len(TMP),nb_regionI], dtype=dty)')
            exec('AREA_' + sim + 'PR = np.zeros([len(TMP),nb_regionI], dtype=dty)')
            for i in range(1, 114 + 1):
                exec(
                    var + '_' + sim + 'PR[:,regionI_index[i]] += TMP[:,i-1]*TMP2[:len(TMP),i-1]')
                exec('AREA_' + sim + 'PR[:,regionI_index[i]] += TMP2[:len(TMP),i-1]')
            exec(var + '_' + sim + 'PR /= AREA_' + sim + 'PR')
            exec(
                var + '_' + sim + 'PR[np.isnan(' + var + '_' + sim + 'PR)|np.isinf(' + var + '_' + sim + 'PR)] = 0')

# precipitation patterns based on 'historical' and 'rcp'
elif (mod_PRECpattern == 'hist&RCPs'):
    for sim in ['ctrl', 'hist', 'rcp26', 'rcp45', 'rcp60', 'rcp85']:
        # global
        if os.path.isfile('data/ClimReg_CMIP5/#DATA.ClimReg_' + mod_PRECresp + '.' + prd[
            sim] + '_(3var).' + sim + '_global.csv'):
            TMP = np.array([line for line in csv.reader(open(
                'data/ClimReg_CMIP5/#DATA.ClimReg_' + mod_PRECresp + '.' + prd[
                    sim] + '_(3var).' + sim + '_global.csv', 'r'))][1:], dtype=dty)
            lgd = [line for line in csv.reader(open(
                'data/ClimReg_CMIP5/#DATA.ClimReg_' + mod_PRECresp + '.' + prd[
                    sim] + '_(3var).' + sim + '_global.csv', 'r'))][0]
            exec('gyp_' + sim + 'PR = TMP[:,lgd.index("pr")]')
        # local
        for var in ['pr']:
            if os.path.isfile(
                    'data/ClimReg_CMIP5/#DATA.ClimReg_' + mod_PRECresp + '.' + prd[
                        sim] + '_114reg1.' + sim + '_' + var + '.csv'):
                TMP = np.array([line for line in csv.reader(open(
                    'data/ClimReg_CMIP5/#DATA.ClimReg_' + mod_PRECresp + '.' + prd[
                        sim] + '_114reg1.' + sim + '_' + var + '.csv', 'r'))], dtype=dty)
                TMP2 = np.array([line for line in csv.reader(open(
                    'data/ClimReg_CMIP5/#DATA.ClimReg_' + mod_PRECresp + '.251yr_114reg1.AREA.csv',
                    'r'))], dtype=dty)
                exec(var + '_' + sim + 'PR = np.zeros([len(TMP),nb_regionI], dtype=dty)')
                exec('AREA_' + sim + 'PR = np.zeros([len(TMP),nb_regionI], dtype=dty)')
                for i in range(1, 114 + 1):
                    exec(
                        var + '_' + sim + 'PR[:,regionI_index[i]] += TMP[:,i-1]*TMP2[:len(TMP),i-1]')
                    exec('AREA_' + sim + 'PR[:,regionI_index[i]] += TMP2[:len(TMP),i-1]')
                exec(var + '_' + sim + 'PR /= AREA_' + sim + 'PR')
                exec(
                    var + '_' + sim + 'PR[np.isnan(' + var + '_' + sim + 'PR)|np.isinf(' + var + '_' + sim + 'PR)] = 0')

# aggregate all experiments
for VAR in ['gyp', 'pr']:
    if (mod_PRECpattern == '4xCO2'):
        exec(VAR + '_allPR = np.array(list(' + VAR + '_quadPR), dtype=dty)')
    elif (mod_PRECpattern == 'hist&RCPs'):
        exec(VAR + '_allPR = []')
        nrcp = 0
        for sim in ['rcp26', 'rcp45', 'rcp60', 'rcp85']:
            if os.path.isfile(
                    'data/ClimReg_CMIP5/#DATA.ClimReg_' + mod_PRECresp + '.' + prd[
                        sim] + '_(3var).' + sim + '_global.csv'):
                nrcp += 1
                exec(
                    VAR + '_allPR += list(' + VAR + '_histPR)+list(' + VAR + '_' + sim + 'PR)')
        exec('test = ' + VAR + '_allPR')
        if (test == []):
            exec(VAR + '_allPR += list(' + VAR + '_histPR)')
        exec(VAR + '_allPR = np.array(' + VAR + '_allPR, dtype=dty)')
# decadal means
for VAR in ['gyp', 'pr']:
    if (mod_PRECpattern == '4xCO2'):
        exec(
            VAR + '_decPR= np.zeros([lng[mod_PRECresp]-10]+list(np.shape(' + VAR + '_allPR)[1:]), dtype=dty)')
        for t in range(lng[mod_PRECresp] - 10):
            exec(VAR + '_decPR[t,...] = np.mean(' + VAR + '_allPR[t:t+10,...],0)')
    elif (mod_PRECpattern == 'hist&RCPs'):
        exec(
            VAR + '_decPR= np.zeros([(156-10)+nrcp*(95-10)]+list(np.shape(' + VAR + '_allPR)[1:]), dtype=dty)')
        for t in range(156 - 10):
            exec(VAR + '_decPR[t,...] = np.mean(' + VAR + '_allPR[t:t+10,...],0)')
        for n in range(nrcp):
            for t in range(95 - 10):
                exec(
                    VAR + '_decPR[(156-10)+n*(95-10)+t,...] = np.mean(' + VAR + '_allPR[(156-10)+n*(95-10)+t:(156-10)+n*(95-10)+t+10,...],0)')

# definition of parameter
# scaling of local yearly precipitations over gst {.}
w_reg_lyp = np.zeros([nb_regionI], dtype=dty)

# fit of parameter
for i in range(nb_regionI):
    diff = pr_decPR[:, i] - np.mean(pr_ctrlPR[:, i], 0)


    def err(var):
        clim = var[0] * (gyp_decPR - np.mean(gyp_ctrlPR, 0))
        return np.sum((diff - clim) ** 2)


    [w_reg_lyp[i]] = fmin(err, [1], disp=False)









################################################################################
# OSCAR FORMAT
################################################################################


##################################################
#   A. FORMAT DRIVERS
##################################################
print 'FORMATING'

# remove attribution axis
# drivers
for VAR in ['EFF','ECH4','EN2O']+['LUC','HARV','SHIFT']+['EHFC','EPFC','EODS']+['ENOX','ECO','EVOC','ESO2','ENH3','EOC','EBC']+['RFcon','RFvolc','RFsolar']:
    exec(VAR+' = np.sum(np.sum(np.sum('+VAR+',3),2),1)')
# parameters
for VAR in ['ECH4','EN2O']+['ENOX','ECO','EVOC','ESO2','ENH3','EOC','EBC']:
    exec(VAR+'_0 = np.sum(np.sum(np.sum('+VAR+'_0,2),1),0)')

################################################################################
# OSCAR FUNCTION
################################################################################


import csv

import numpy as np
import matplotlib.pyplot as plt

from matplotlib.font_manager import FontProperties


##################################################
#   1. OSCAR LITE
##################################################

def OSCAR_lite(p=p, fT=fT, \
               EFF=EFF, ECH4=ECH4, EN2O=EN2O, \
               LUC=LUC, HARV=HARV, SHIFT=SHIFT, \
               EHFC=EHFC, EPFC=EPFC, EODS=EODS, \
               ENOX=ENOX, ECO=ECO, EVOC=EVOC, ESO2=ESO2, ENH3=ENH3, EOC=EOC, EBC=EBC, \
               RFcon=RFcon, RFvolc=RFvolc, RFsolar=RFsolar,
               force_CO2=False, force_GHG=False, force_halo=False, force_RF=False,
               force_RFs=False, force_clim=False, \
               var_output=['ELUC', 'OSNK', 'LSNK', 'D_CO2', 'RF', 'D_gst'], \
               plot=[]):
    # ===============
    # A. DEFINITIONS
    # ===============

    # plot variables
    var_plot = []
    if plot == 'all' or plot == 'CO2' or 'CO2' in plot:
        var_plot += ['D_CO2', 'OSNK', 'LSNK', 'ELUC', 'D_AREA', 'D_npp', 'D_efire',
                     'D_fmort', 'D_rh1', 'D_fmet', 'D_rh2', 'D_FIN', 'D_FOUT', 'D_FCIRC',
                     'EFIRE_luc', 'FMORT_luc', 'RH1_luc', 'FMET_luc', 'RH2_luc',
                     'EHWP1_luc', 'EHWP2_luc', 'EHWP3_luc']
    if plot == 'all' or plot == 'CH4' or 'CH4' in plot:
        var_plot += ['D_CH4', 'D_OHSNK_CH4', 'D_HVSNK_CH4', 'D_XSNK_CH4', 'D_EWET',
                     'D_EBB_CH4']
    if plot == 'all' or plot == 'N2O' or 'N2O' in plot:
        var_plot += ['D_N2O', 'D_HVSNK_N2O', 'D_EBB_N2O']
    if plot == 'all' or plot == 'O3' or 'O3' in plot:
        var_plot += ['D_O3t', 'D_O3s', 'D_EESC', 'D_N2O_lag', 'D_gst']
    if plot == 'all' or plot == 'AER' or 'AER' in plot:
        var_plot += ['D_SO4', 'D_POA', 'D_BC', 'D_NO3', 'D_SOA', 'D_AERh', 'RF_SO4',
                     'RF_POA', 'RF_BC', 'RF_NO3', 'RF_SOA', 'RF_cloud']
    if plot == 'all' or plot == 'clim' or 'clim' in plot:
        var_plot += ['RF', 'D_gst', 'D_gyp', 'RF_CO2', 'RF_CH4', 'RF_H2Os', 'RF_N2O',
                     'RF_halo', 'RF_O3t', 'RF_O3s', 'RF_SO4', 'RF_POA', 'RF_BC', 'RF_NO3',
                     'RF_SOA', 'RF_cloud', 'RF_BCsnow', 'RF_LCC']

    # save variables
    var_timeseries = list(set(var_output) | set(var_plot))
    for var in var_timeseries:
        # global variables
        if var in ['D_mld', 'D_dic', 'D_pH'] \
                or var in ['OSNK', 'LSNK', 'D_FOXI_CH4', 'D_OHSNK_CH4', 'D_HVSNK_CH4',
                           'D_XSNK_CH4', 'D_HVSNK_N2O'] \
                or var in ['D_kOH', 'D_hv'] \
                or var in ['D_O3t', 'D_EESC', 'D_O3s', 'D_SO4', 'D_POA', 'D_BC', 'D_NO3',
                           'D_SOA', 'D_AERh'] \
                or var in ['D_CO2', 'D_CH4', 'D_CH4_lag', 'D_N2O', 'D_N2O_lag'] \
                or var in ['RF', 'RF_warm', 'RF_atm', 'RF_CO2', 'RF_CH4', 'RF_H2Os',
                           'RF_N2O', 'RF_halo', 'RF_O3t', 'RF_O3s', 'RF_SO4', 'RF_POA',
                           'RF_BC', 'RF_NO3', 'RF_SOA', 'RF_cloud', 'RF_BCsnow', 'RF_LCC'] \
                or var in ['D_gst', 'D_sst', 'D_gyp', 'D_OHC']:
            exec(var + '_t = np.zeros([ind_final+1],dtype=dty)')
        # (region) variables
        if var in ['ELUC', 'D_AWET', 'D_EWET', 'D_ewet', 'D_EBB_CO2', 'D_EBB_CH4',
                   'D_EBB_N2O', 'D_EBB_NOX', 'D_EBB_CO', 'D_EBB_VOC', 'D_EBB_SO2',
                   'D_EBB_NH3', 'D_EBB_OC', 'D_EBB_BC', 'D_lst', 'D_lyp']:
            exec(var + '_t = np.zeros([ind_final+1,nb_regionI],dtype=dty)')
            # (region)*(biome) variables
        if var in ['D_AREA', 'D_npp', 'D_efire', 'D_fmort', 'D_rh1', 'D_fmet', 'D_rh2',
                   'D_cveg', 'D_csoil1', 'D_csoil2']:
            exec(var + '_t = np.zeros([ind_final+1,nb_regionI,nb_biome],dtype=dty)')
        # (region)*(biome)*(biome)*(age) variables
        if var in ['EFIRE_luc', 'FMORT_luc', 'RH1_luc', 'FMET_luc', 'RH2_luc',
                   'EHWP1_luc', 'EHWP2_luc', 'EHWP3_luc'] \
                or var in ['CVEG_luc', 'CSOIL1_luc', 'CSOIL2_luc', 'CHWP1_luc',
                           'CHWP2_luc', 'CHWP3_luc']:
            exec(
                var + '_t = np.zeros([ind_final+1,nb_regionI,nb_biome,nb_biome,ind_final+1],dtype=dty)')
            # (obox) variables
        if var in ['D_FIN', 'D_FOUT', 'D_FCIRC', 'D_CSURF']:
            exec(var + '_t = np.zeros([ind_final+1,nb_obox],dtype=dty)')
        # (species) variables
        if var in ['D_HFC', 'D_HFC_lag', 'D_OHSNK_HFC', 'D_HVSNK_HFC', 'D_XSNK_HFC']:
            exec(var + '_t = np.zeros([ind_final+1,nb_HFC],dtype=dty)')
        if var in ['D_PFC', 'D_PFC_lag', 'D_OHSNK_PFC', 'D_HVSNK_PFC', 'D_XSNK_PFC']:
            exec(var + '_t = np.zeros([ind_final+1,nb_PFC],dtype=dty)')
        if var in ['D_ODS', 'D_ODS_lag', 'D_OHSNK_ODS', 'D_HVSNK_ODS', 'D_XSNK_ODS']:
            exec(var + '_t = np.zeros([ind_final+1,nb_ODS],dtype=dty)')
        # (regionPF) variables
        if var in ['pthaw', 'pthaw_bar'] \
                or var in ['FTHAW', 'ETHAW1', 'ETHAW2', 'ETHAW3', 'EPF_CO2', 'EPF_CH4',
                           'EPF'] \
                or var in ['CTHAW1', 'CTHAW2', 'CTHAW3', 'D_CFROZ']:
            exec(var + '_t = np.zeros([ind_final+1,nb_regionPF], dtype=dty)')

    # run variables
    # ocean
    D_dic = np.array([0], dtype=dty)
    D_CSURF = np.zeros([nb_obox], dtype=dty)
    # land
    for var in ['D_AREA', 'D_cveg', 'D_csoil1', 'D_csoil2']:
        exec(var + ' = np.zeros([nb_regionI,nb_biome],dtype=dty)')
    # land-use
    for var in ['CVEG_luc', 'CSOIL1_luc', 'CSOIL2_luc', 'CHWP1_luc', 'CHWP2_luc',
                'CHWP3_luc']:
        exec(var + ' = np.zeros([nb_regionI,nb_biome,nb_biome,ind_final+1],dtype=dty)')
    # atmosphere
    for var in ['D_CO2', 'D_CH4', 'D_CH4_lag', 'D_N2O', 'D_N2O_lag', 'D_EESC', 'D_O3s']:
        exec(var + ' = np.array([0],dtype=dty)')
    for var in ['D_HFC', 'D_HFC_lag', 'D_PFC', 'D_PFC_lag', 'D_ODS', 'D_ODS_lag']:
        exec(var + ' = np.zeros([nb_' + var[2:2 + 3] + '],dtype=dty)')
    # climate
    for var in ['D_gst', 'D_gst0', 'D_sst', 'D_gyp', 'D_OHC']:
        exec(var + ' = np.array([0],dtype=dty)')
    for var in ['D_lst', 'D_lyp']:
        exec(var + ' = np.zeros([nb_regionI],dtype=dty)')
    # permafrost
    for var in ['pthaw', 'CTHAW1', 'CTHAW2', 'CTHAW3', 'D_CFROZ']:
        exec(var + ' = np.zeros([nb_regionPF], dtype=dty)')

    # =======
    # B. RUN
    # =======

    for t in range(1, ind_final + 1):
        for tt in range(p):

            # ---------
            # 1. OCEAN
            # ---------

            # structure
            D_mld = mld_0 * alpha_mld * (np.exp(gamma_mld * fT * D_sst) - 1)
            # fluxes
            D_FIN = p_circ * v_fg * alpha_CO2 * D_CO2
            D_FOUT = p_circ * v_fg * alpha_CO2 * f_pCO2(D_dic, fT * D_sst)
            D_FCIRC = D_CSURF * (1 / tau_circ)
            OSNK = np.sum(D_FOUT - D_FIN)
            # stocks
            D_CSURF += (p ** -1) * (D_FIN - D_FOUT - D_FCIRC)
            D_dic = alpha_dic * np.sum(D_CSURF) / (1 + D_mld / mld_0)

            # --------
            # 2. LAND
            # --------

            # land-cover
            D_AREA += (p ** -1) * (np.sum(LUC[t], 1) - np.sum(LUC[t], 2))
            D_AWET = AWET_0 * (
                        gamma_wetT * fT * D_lst + gamma_wetP * fT * D_lyp + gamma_wetC * fT * D_CO2)
            # factors
            D_k_igni = gamma_igniT * fT * D_lst[:, np.newaxis] + gamma_igniP * fT * D_lyp[
                                                                                    :,
                                                                                    np.newaxis] + gamma_igniC * fT * D_CO2
            D_k_rho = f_rho(fT * D_lst[:, np.newaxis], fT * D_lyp[:, np.newaxis])
            # fluxes
            D_npp = npp_0 * f_npp(D_CO2, fT * D_lst[:, np.newaxis],
                                  fT * D_lyp[:, np.newaxis])
            D_efire = igni_0 * ((1 + D_k_igni) * (cveg_0 + D_cveg) - cveg_0)
            D_fmort = mu_0 * D_cveg
            D_rh1 = rho1_0 * ((1 + D_k_rho) * (csoil1_0 + D_csoil1) - csoil1_0)
            D_fmet = k_met * D_rh1
            D_rh2 = rho2_0 * ((1 + D_k_rho) * (csoil2_0 + D_csoil2) - csoil2_0)
            D_ewet = ewet_0 * np.nan_to_num(
                np.sum(p_wet * D_csoil1, 1) / np.sum(p_wet * csoil1_0, 1))
            LSNK = np.sum((D_rh1 + D_rh2 + D_efire - D_npp) * (AREA_0 + D_AREA))
            D_EWET = ewet_0 * D_AWET + D_ewet * AWET_0 + D_ewet * D_AWET
            # stocks
            D_cveg += (p ** -1) * (D_npp - D_fmort - D_efire)
            D_csoil1 += (p ** -1) * (D_fmort - D_fmet - D_rh1)
            D_csoil2 += (p ** -1) * (D_fmet - D_rh2)

            # ------------
            # 3. LAND-USE
            # ------------

            # initialization
            # land-use change
            for b1 in range(nb_biome):
                for b2 in range(nb_biome):
                    CVEG_luc[:, b1, b2, t] += (p ** -1) * -(cveg_0 + D_cveg)[:, b2] * LUC[
                                                                                      t,
                                                                                      :,
                                                                                      b1,
                                                                                      b2]
                    CSOIL1_luc[:, b1, b2, t] += (p ** -1) * (
                                (csoil1_0 + D_csoil1)[:, b1] - (csoil1_0 + D_csoil1)[:,
                                                               b2]) * LUC[t, :, b1, b2]
                    CSOIL2_luc[:, b1, b2, t] += (p ** -1) * (
                                (csoil2_0 + D_csoil2)[:, b1] - (csoil2_0 + D_csoil2)[:,
                                                               b2]) * LUC[t, :, b1, b2]
                    CSOIL1_luc[:, b1, b2, t] += (p ** -1) * (cveg_0 + D_cveg)[:, b1] * (
                                p_AGB[:, b1] * p_HWP0[:, b1] + (1 - p_AGB[:, b1])) * LUC[
                                                                                     t, :,
                                                                                     b1,
                                                                                     b2]
                    CHWP1_luc[:, b1, b2, t] += (p ** -1) * (cveg_0 + D_cveg)[:,
                                                           b1] * p_AGB[:, b1] * p_HWP1[:,
                                                                                b1] * LUC[
                                                                                      t,
                                                                                      :,
                                                                                      b1,
                                                                                      b2]
                    CHWP2_luc[:, b1, b2, t] += (p ** -1) * (cveg_0 + D_cveg)[:,
                                                           b1] * p_AGB[:, b1] * p_HWP2[:,
                                                                                b1] * LUC[
                                                                                      t,
                                                                                      :,
                                                                                      b1,
                                                                                      b2]
                    CHWP3_luc[:, b1, b2, t] += (p ** -1) * (cveg_0 + D_cveg)[:,
                                                           b1] * p_AGB[:, b1] * p_HWP3[:,
                                                                                b1] * LUC[
                                                                                      t,
                                                                                      :,
                                                                                      b1,
                                                                                      b2]
            # harvest
            for b in range(nb_biome):
                CVEG_luc[:, b, b, t] += (p ** -1) * -HARV[t, :, b]
                CSOIL1_luc[:, b, b, t] += (p ** -1) * p_HWP0[:, b] * HARV[t, :, b]
                CHWP1_luc[:, b, b, t] += (p ** -1) * p_HWP1[:, b] * HARV[t, :, b]
                CHWP2_luc[:, b, b, t] += (p ** -1) * p_HWP2[:, b] * HARV[t, :, b]
                CHWP3_luc[:, b, b, t] += (p ** -1) * p_HWP3[:, b] * HARV[t, :, b]
            # shifting cultivation
            for b1 in range(nb_biome):
                for b2 in range(b1, nb_biome):
                    CVEG_luc[:, b1, b2, t] += (p ** -1) * -(cveg_0 + D_cveg)[:, b2] * (
                                1 - np.exp(-mu_0[:, b2] * tau_shift)) * SHIFT[t, :, b1,
                                                                        b2]
                    CSOIL1_luc[:, b1, b2, t] += (p ** -1) * (cveg_0 + D_cveg)[:, b1] * (
                                1 - np.exp(-mu_0[:, b1] * tau_shift)) * (
                                                            p_AGB[:, b1] * p_HWP0[:,
                                                                           b1] + (
                                                                        1 - p_AGB[:,
                                                                            b1])) * SHIFT[
                                                                                    t, :,
                                                                                    b1,
                                                                                    b2]
                    CHWP1_luc[:, b1, b2, t] += (p ** -1) * (cveg_0 + D_cveg)[:, b1] * (
                                1 - np.exp(-mu_0[:, b1] * tau_shift)) * p_AGB[:,
                                                                        b1] * p_HWP1[:,
                                                                              b1] * SHIFT[
                                                                                    t, :,
                                                                                    b1,
                                                                                    b2]
                    CHWP2_luc[:, b1, b2, t] += (p ** -1) * (cveg_0 + D_cveg)[:, b1] * (
                                1 - np.exp(-mu_0[:, b1] * tau_shift)) * p_AGB[:,
                                                                        b1] * p_HWP2[:,
                                                                              b1] * SHIFT[
                                                                                    t, :,
                                                                                    b1,
                                                                                    b2]
                    CHWP3_luc[:, b1, b2, t] += (p ** -1) * (cveg_0 + D_cveg)[:, b1] * (
                                1 - np.exp(-mu_0[:, b1] * tau_shift)) * p_AGB[:,
                                                                        b1] * p_HWP3[:,
                                                                              b1] * SHIFT[
                                                                                    t, :,
                                                                                    b1,
                                                                                    b2]
                    CVEG_luc[:, b2, b1, t] += (p ** -1) * -(cveg_0 + D_cveg)[:, b1] * (
                                1 - np.exp(-mu_0[:, b1] * tau_shift)) * SHIFT[t, :, b1,
                                                                        b2]
                    CSOIL1_luc[:, b2, b1, t] += (p ** -1) * (cveg_0 + D_cveg)[:, b2] * (
                                1 - np.exp(-mu_0[:, b2] * tau_shift)) * (
                                                            p_AGB[:, b2] * p_HWP0[:,
                                                                           b2] + (
                                                                        1 - p_AGB[:,
                                                                            b2])) * SHIFT[
                                                                                    t, :,
                                                                                    b1,
                                                                                    b2]
                    CHWP1_luc[:, b2, b1, t] += (p ** -1) * (cveg_0 + D_cveg)[:, b2] * (
                                1 - np.exp(-mu_0[:, b2] * tau_shift)) * p_AGB[:,
                                                                        b2] * p_HWP1[:,
                                                                              b2] * SHIFT[
                                                                                    t, :,
                                                                                    b1,
                                                                                    b2]
                    CHWP2_luc[:, b2, b1, t] += (p ** -1) * (cveg_0 + D_cveg)[:, b2] * (
                                1 - np.exp(-mu_0[:, b2] * tau_shift)) * p_AGB[:,
                                                                        b2] * p_HWP2[:,
                                                                              b2] * SHIFT[
                                                                                    t, :,
                                                                                    b1,
                                                                                    b2]
                    CHWP3_luc[:, b2, b1, t] += (p ** -1) * (cveg_0 + D_cveg)[:, b2] * (
                                1 - np.exp(-mu_0[:, b2] * tau_shift)) * p_AGB[:,
                                                                        b2] * p_HWP3[:,
                                                                              b2] * SHIFT[
                                                                                    t, :,
                                                                                    b1,
                                                                                    b2]

            # fluxes
            # book-keeping model
            NPP_luc = 0 * CVEG_luc
            EFIRE_luc = (igni_0 * (1 + D_k_igni))[:, np.newaxis, :, np.newaxis] * CVEG_luc
            FMORT_luc = mu_0[:, np.newaxis, :, np.newaxis] * CVEG_luc
            RH1_luc = (rho1_0 * (1 + D_k_rho))[:, np.newaxis, :, np.newaxis] * CSOIL1_luc
            FMET_luc = k_met * RH1_luc
            RH2_luc = (rho2_0 * (1 + D_k_rho))[:, np.newaxis, :, np.newaxis] * CSOIL2_luc
            EHWP1_luc = np.zeros([nb_regionI, nb_biome, nb_biome, ind_final + 1],
                                 dtype=dty)
            EHWP1_luc[:, :, :, :t + 1] = r_HWP1[np.newaxis, np.newaxis, np.newaxis,
                                         t::-1] * CHWP1_luc[:, :, :, :t + 1]
            EHWP2_luc = np.zeros([nb_regionI, nb_biome, nb_biome, ind_final + 1],
                                 dtype=dty)
            EHWP2_luc[:, :, :, :t + 1] = r_HWP2[np.newaxis, np.newaxis, np.newaxis,
                                         t::-1] * CHWP2_luc[:, :, :, :t + 1]
            EHWP3_luc = np.zeros([nb_regionI, nb_biome, nb_biome, ind_final + 1],
                                 dtype=dty)
            EHWP3_luc[:, :, :, :t + 1] = r_HWP3[np.newaxis, np.newaxis, np.newaxis,
                                         t::-1] * CHWP3_luc[:, :, :, :t + 1]
            ELUC = np.sum(np.sum(np.sum(
                RH1_luc + RH2_luc + EFIRE_luc + EHWP1_luc + EHWP2_luc + EHWP3_luc - NPP_luc,
                3), 2), 1)
            # biomass burning
            for VAR in ['CO2', 'CH4', 'N2O', 'NOX', 'CO', 'VOC', 'SO2', 'NH3', 'OC',
                        'BC']:
                exec(
                    'D_EBB_' + VAR + ' = np.sum( alpha_BB_' + VAR + '*(igni_0*cveg_0*D_AREA + D_efire*AREA_0 + D_efire*D_AREA) ,1)')
                exec(
                    'D_EBB_' + VAR + ' += p_HWP1_bb * np.sum(np.sum(np.sum( alpha_BB_' + VAR + '[:,:,np.newaxis,np.newaxis]*EHWP1_luc ,3),2),1)')
                exec(
                    'D_EBB_' + VAR + ' += np.sum(np.sum(np.sum( alpha_BB_' + VAR + '[:,np.newaxis,:,np.newaxis]*EFIRE_luc ,3),2),1)')

            # stocks
            CVEG_luc += (p ** -1) * (NPP_luc - FMORT_luc - EFIRE_luc)
            CSOIL1_luc += (p ** -1) * (FMORT_luc - FMET_luc - RH1_luc)
            CSOIL2_luc += (p ** -1) * (FMET_luc - RH2_luc)
            CHWP1_luc += (p ** -1) * -EHWP1_luc
            CHWP2_luc += (p ** -1) * -EHWP2_luc
            CHWP3_luc += (p ** -1) * -EHWP3_luc

            # ---------------
            # 3b. PERMAFROST
            # ---------------

            # factors
            rD_rhoPF = np.exp(
                w_rhoPF * gamma_rhoPF1 * w_reg_lstPF * fT * D_gst + w_rhoPF * gamma_rhoPF2 * (
                        w_reg_lstPF * fT * D_gst) ** 2) - 1
            # fraction thawed
            pthaw_bar = -pthaw_min + (1 + pthaw_min) / (
                        1 + ((1 / pthaw_min + 1) ** k_pthaw - 1) * np.exp(
                    -gamma_pthaw * k_pthaw * w_reg_lstPF * fT * D_gst)) ** (1 / k_pthaw)
            d_pthaw = f_v_PF(pthaw_bar, pthaw) * (pthaw_bar - pthaw)
            pthaw += (p ** -1) * d_pthaw
            # fluxes
            FTHAW = CFROZ_0 * d_pthaw
            ETHAW1 = 1 / tau_PF1 * (1 + rD_rhoPF) * CTHAW1
            ETHAW2 = 1 / tau_PF2 * (1 + rD_rhoPF) * CTHAW2
            ETHAW3 = 1 / tau_PF3 * (1 + rD_rhoPF) * CTHAW3
            EPF_CO2 = (1 - p_PF_CH4) * (ETHAW1 + ETHAW2 + ETHAW3 + p_PF_inst * FTHAW)
            EPF_CH4 = 1000. * p_PF_CH4 * (ETHAW1 + ETHAW2 + ETHAW3 + p_PF_inst * FTHAW)
            EPF = EPF_CO2 + 0.001 * EPF_CH4
            # stocks
            D_CFROZ -= (p ** -1) * FTHAW
            CTHAW1 += (p ** -1) * (p_PF1 * (1 - p_PF_inst) * FTHAW - ETHAW1)
            CTHAW2 += (p ** -1) * (p_PF2 * (1 - p_PF_inst) * FTHAW - ETHAW2)
            CTHAW3 += (p ** -1) * (p_PF3 * (1 - p_PF_inst) * FTHAW - ETHAW3)

            # -------------
            # 4. CHEMISTRY
            # -------------

            # factors
            D_kOH = f_kOH(D_CH4, D_O3s, fT * D_gst, np.sum(ENOX[t] + D_EBB_NOX),
                          np.sum(ECO[t] + D_EBB_CO), np.sum(EVOC[t] + D_EBB_VOC))
            D_hv = f_hv(D_N2O_lag, D_EESC, fT * D_gst)
            # fluxes
            D_OHSNK_CH4 = -alpha_CH4 / tau_CH4_OH * (
                        CH4_0 * D_kOH + D_CH4 + D_kOH * D_CH4)
            D_HVSNK_CH4 = -alpha_CH4 / tau_CH4_hv * (
                        CH4_0 * D_hv + D_CH4_lag + D_hv * D_CH4_lag)
            D_XSNK_CH4 = -alpha_CH4 * (1 / tau_CH4_soil + 1 / tau_CH4_ocean) * D_CH4
            D_FOXI_CH4 = -0.001 * (1.0 * np.sum(ECH4[t]) + np.sum(D_EBB_CH4) + np.sum(
                D_EWET) + D_OHSNK_CH4 + D_HVSNK_CH4 + D_XSNK_CH4)
            D_HVSNK_N2O = -alpha_N2O / tau_N2O_hv * (
                        N2O_0 * D_hv + D_N2O_lag + D_hv * D_N2O_lag)
            for VAR in ['HFC', 'PFC', 'ODS']:
                exec(
                    'D_OHSNK_' + VAR + ' = -alpha_' + VAR + '/tau_' + VAR + '_OH * (' + VAR + '_0*D_kOH + D_' + VAR + ' + D_kOH*D_' + VAR + ')')
                exec(
                    'D_HVSNK_' + VAR + ' = -alpha_' + VAR + '/tau_' + VAR + '_hv * (' + VAR + '_0*D_hv + D_' + VAR + '_lag + D_hv*D_' + VAR + '_lag)')
                exec(
                    'D_XSNK_' + VAR + ' = -alpha_' + VAR + '/tau_' + VAR + '_othr * D_' + VAR)
            # stocks
            D_O3t = chi_O3t_CH4 * np.log(1 + D_CH4 / CH4_0) + Gamma_O3t * fT * D_gst
            D_O3t += chi_O3t_NOX * np.sum(
                w_reg_NOX * np.sum(p_reg4 * (ENOX[t] + D_EBB_NOX)[:, np.newaxis], 0))
            D_O3t += chi_O3t_CO * np.sum(
                w_reg_CO * np.sum(p_reg4 * (ECO[t] + D_EBB_CO)[:, np.newaxis], 0))
            D_O3t += chi_O3t_VOC * np.sum(
                w_reg_VOC * np.sum(p_reg4 * (EVOC[t] + D_EBB_VOC)[:, np.newaxis], 0))
            D_EESC = np.sum(f_fracrel(tau_lag) * (n_Cl + alpha_Br * n_Br) * D_ODS_lag)
            D_O3s = chi_O3s_EESC * D_EESC + chi_O3s_N2O * D_N2O_lag * (
                        1 - D_EESC / EESC_x) + Gamma_O3s * fT * D_gst
            D_SO4 = alpha_SO4 * tau_SO2 * np.sum(
                w_reg_SO2 * np.sum(p_reg4 * (ESO2[t] + D_EBB_SO2)[:, np.newaxis],
                                   0)) + alpha_SO4 * tau_DMS * 0 + Gamma_SO4 * fT * D_gst
            D_POA = tau_OMff * alpha_POM * np.sum(
                w_reg_OC * np.sum(p_reg4 * (EOC[t])[:, np.newaxis],
                                  0)) + tau_OMbb * alpha_POM * np.sum(
                D_EBB_OC) + Gamma_POA * fT * D_gst
            D_BC = tau_BCff * np.sum(w_reg_BC * np.sum(p_reg4 * (EBC[t])[:, np.newaxis],
                                                       0)) + tau_BCbb * np.sum(
                D_EBB_BC) + Gamma_BC * fT * D_gst
            D_NO3 = alpha_NO3 * tau_NOX * np.sum(
                ENOX[t] + D_EBB_NOX) + alpha_NO3 * tau_NH3 * np.sum(
                ENH3[t] + D_EBB_NH3) + Gamma_NO3 * fT * D_gst
            D_SOA = tau_VOC * np.sum(
                EVOC[t] + D_EBB_VOC) + tau_BVOC * 0 + Gamma_SOA * fT * D_gst
            D_DUST = 0 * (tau_DUST * 0 + Gamma_DUST * fT * D_gst)
            D_SALT = 0 * (tau_SALT * 0 + Gamma_SALT * fT * D_gst)
            D_AERh = solub_SO4 * D_SO4 + solub_POA * D_POA + solub_BC * D_BC + solub_NO3 * D_NO3 + solub_SOA * D_SOA + solub_DUST * D_DUST + solub_SALT * D_SALT

            # --------------
            # 5. ATMOSPHERE
            # --------------

            # stocks
            D_CO2 += (p ** -1) * (1 / alpha_CO2) * (
                        np.sum(EFF[t]) + np.sum(ELUC) + LSNK + OSNK + D_FOXI_CH4 + np.sum(
                    EPF_CO2))
            D_CH4 += (p ** -1) * (1 / alpha_CH4) * (
                        np.sum(ECH4[t]) + np.sum(D_EBB_CH4) + np.sum(D_EWET) + np.sum(
                    EPF_CH4) + D_OHSNK_CH4 + D_HVSNK_CH4 + D_XSNK_CH4)
            D_N2O += (p ** -1) * (1 / alpha_N2O) * (
                        np.sum(EN2O[t]) + np.sum(D_EBB_N2O) + D_HVSNK_N2O)
            D_HFC += (p ** -1) * (1 / alpha_HFC) * (
                        np.sum(EHFC[t], 0) + D_OHSNK_HFC + D_HVSNK_HFC + D_XSNK_HFC)
            D_PFC += (p ** -1) * (1 / alpha_PFC) * (
                        np.sum(EPFC[t], 0) + D_OHSNK_PFC + D_HVSNK_PFC + D_XSNK_PFC)
            D_ODS += (p ** -1) * (1 / alpha_ODS) * (
                        np.sum(EODS[t], 0) + D_OHSNK_ODS + D_HVSNK_ODS + D_XSNK_ODS)
            for VAR in ['CH4', 'N2O', 'HFC', 'PFC', 'ODS']:
                exec(
                    'D_' + VAR + '_lag += (p**-1) * ((1/tau_lag)*D_' + VAR + ' - (1/tau_lag)*D_' + VAR + '_lag)')

            # FORCE
            if force_CO2:
                D_CO2 = D_CO2_force[t]

            if force_GHG:
                D_CO2 = D_CO2_force[t]
                D_CH4 = D_CH4_force[t]
                D_N2O = D_N2O_force[t]

            if force_halo:
                D_HFC[:] = D_HFC_force[t]
                D_PFC[:] = D_PFC_force[t]
                D_ODS[:] = D_ODS_force[t]

            # -----------
            # 6. CLIMATE
            # -----------

            # fluxes
            # per component
            RF_CO2 = f_RF_CO2(D_CO2)
            RF_CH4 = f_RF_CH4(D_CH4) - (
                        f_RF_overlap(D_CH4, D_N2O) - f_RF_overlap(0, D_N2O))
            RF_H2Os = f_RF_H2Os(D_CH4_lag)
            RF_N2O = f_RF_N2O(D_N2O) - (
                        f_RF_overlap(D_CH4, D_N2O) - f_RF_overlap(D_CH4, 0))
            RF_halo = np.sum(radeff_HFC * D_HFC) + np.sum(radeff_PFC * D_PFC) + np.sum(
                radeff_ODS * D_ODS)
            for VAR in ['O3t', 'O3s', 'SO4', 'POA', 'BC', 'NO3', 'SOA', 'DUST', 'SALT']:
                exec('RF_' + VAR + ' = radeff_' + VAR + '*D_' + VAR)
            RF_cloud = k_BC_adjust * RF_BC + Phi_0 * np.log(1 + D_AERh / AERh_0)
            RF_BCsnow = radeff_BCsnow * np.sum(
                w_reg_BCsnow * np.sum(p_reg9 * (EBC[t] + D_EBB_BC)[:, np.newaxis], 0))
            RF_LCC = np.sum(alpha_LCC * D_AREA)

            # FORCE
            if force_RFs:
                for VAR in ['CO2', 'CH4', 'H2Os', 'N2O', 'halo'] + ['O3t', 'O3s', 'SO4',
                                                                    'POA', 'BC', 'NO3',
                                                                    'SOA', 'DUST',
                                                                    'SALT'] + ['cloud',
                                                                               'BCsnow',
                                                                               'LCC']:
                    exec('RF_' + VAR + ' = RF_' + VAR + '_force[t]')

            # totals
            RF = RF_CO2 + RF_CH4 + RF_H2Os + RF_N2O + RF_halo + RF_O3t + RF_O3s + RF_SO4 + RF_POA + RF_BC + RF_NO3 + RF_SOA + RF_DUST + RF_SALT + RF_cloud + RF_BCsnow + RF_LCC + \
                 RFcon[t] + RFvolc[t] + RFsolar[t]
            RF_warm = RF_CO2 + RF_CH4 + RF_H2Os + RF_N2O + RF_halo + RF_O3t + RF_O3s + RF_SO4 + RF_POA + RF_BC + RF_NO3 + RF_SOA + RF_DUST + RF_SALT + RF_cloud + warmeff_BCsnow * RF_BCsnow + warmeff_LCC * RF_LCC + \
                      RFcon[t] + warmeff_volc * RFvolc[t] + RFsolar[t]
            RF_atm = p_atm_CO2 * RF_CO2 + p_atm_noCO2 * (
                        RF_CH4 + RF_N2O + RF_halo) + p_atm_O3t * RF_O3t + p_atm_strat * (
                                 RF_O3s + RF_H2Os) + p_atm_scatter * (
                                 RF_SO4 + RF_POA + RF_NO3 + RF_SOA + RF_DUST + RF_SALT +
                                 RFvolc[t]) + p_atm_absorb * RF_BC + p_atm_cloud * (
                                 RF_cloud + RFcon[t]) + p_atm_alb * (
                                 RF_BCsnow + RF_LCC) + p_atm_solar * RFsolar[t]

            # FORCE
            if force_RF:
                RF_warm = RF_force[t] * (RF_warm / RF)
                RF_atm = RF_force[t] * (RF_atm / RF)
                RF = RF_force[t]

            # stocks
            # temperatures
            D_gst += (p ** -1) * (1 / tau_gst) * (
                        lambda_0 * RF_warm - D_gst - theta_0 * (D_gst - D_gst0))
            D_gst0 += (p ** -1) * (1 / tau_gst0) * theta_0 * (D_gst - D_gst0)
            D_sst = w_reg_sst * D_gst
            D_lst = w_reg_lst * D_gst
            # precipitations
            D_gyp = alpha_gyp * D_gst + beta_gyp * RF_atm
            D_lyp = w_reg_lyp * D_gyp
            # ocean
            D_OHC += (p ** -1) * p_OHC * alpha_OHC * (RF - D_gst / lambda_0)
            D_pH = f_pH(D_CO2)

            # FORCE
            if force_clim:
                D_gst = D_gst_force[t]
                D_sst = D_sst_force[t]
                D_lst = D_lst_force[t]
                D_lyp = D_lyp_force[t]

            # -----------
            # Y. SAVE
            # -----------

            for var in var_timeseries:
                exec(var + '_t[t] += (p**-1) * ' + var)

            # ---------
            # Z. TESTS
            # ---------

            if np.isnan(np.sum(D_CO2)):
                print
                'D_CO2 = NaN at t = ' + str(t) + ' and tt = ' + str(tt)
                print
                'OSNK = ' + str(np.sum(OSNK))
                print
                'LSNK = ' + str(np.sum(LSNK))
                print
                'ELUC = ' + str(np.sum(ELUC))
                break
            if np.isnan(np.sum(D_CH4)):
                print
                'D_CH4 = NaN at t = ' + str(t) + ' and tt = ' + str(tt)
                print
                'D_EWET = ' + str(np.sum(D_EWET))
                print
                'D_OHSNK = ' + str(np.sum(D_OHSNK_CH4))
                print
                'D_HVSNK = ' + str(np.sum(D_HVSNK_CH4))
                break
            if np.isnan(np.sum(D_gst)):
                print
                'D_gst = NaN at t = ' + str(t) + ' and tt = ' + str(tt)
                print
                'RF_CO2 = ' + str(np.sum(RF_CO2))
                print
                'RF_CH4 = ' + str(np.sum(RF_CH4))
                print
                'RF_H2Os = ' + str(np.sum(RF_H2Os))
                print
                'RF_N2O = ' + str(np.sum(RF_N2O))
                print
                'RF_halo = ' + str(np.sum(RF_halo))
                print
                'RF_O3t = ' + str(np.sum(RF_O3t))
                print
                'RF_O3s = ' + str(np.sum(RF_O3s))
                print
                'RF_SO4 = ' + str(np.sum(RF_SO4))
                print
                'RF_POA = ' + str(np.sum(RF_POA))
                print
                'RF_BC = ' + str(np.sum(RF_BC))
                print
                'RF_NO3 = ' + str(np.sum(RF_NO3))
                print
                'RF_SOA = ' + str(np.sum(RF_SOA))
                print
                'RF_DUST = ' + str(np.sum(RF_DUST))
                print
                'RF_SALT = ' + str(np.sum(RF_SALT))
                print
                'RF_cloud = ' + str(np.sum(RF_cloud))
                print
                'RF_BCsnow = ' + str(np.sum(RF_BCsnow))
                print
                'RF_LCC = ' + str(np.sum(RF_LCC))
                break

        if np.isnan(np.sum(D_CO2)) | np.isnan(np.sum(D_CH4)) | np.isnan(np.sum(D_gst)):
            for var in var_timeseries:
                if (t < ind_final):
                    exec(var + '_t[t+1:] = np.nan')
            break

    # ===========
    # C. FIGURES
    # ===========

    print
    'PLOTTING', plot
    if plot == 'all' or plot == 'CO2' or 'CO2' in plot:
        plot_CO2(D_CO2_t, OSNK_t, LSNK_t, ELUC_t, EFF, D_AREA_t, D_npp_t, D_efire_t,
                 D_fmort_t, D_rh1_t, D_fmet_t, D_rh2_t, D_FIN_t, D_FOUT_t, D_FCIRC_t,
                 EFIRE_luc_t, FMORT_luc_t, RH1_luc_t, RH2_luc_t, EHWP1_luc_t, EHWP2_luc_t,
                 EHWP3_luc_t)
        plt.savefig('results/plot-CO2.svg')
    if plot == 'all' or plot == 'CH4' or 'CH4' in plot:
        plot_CH4(D_CH4_t, D_OHSNK_CH4_t, D_HVSNK_CH4_t, D_XSNK_CH4_t, D_EWET_t,
                 D_EBB_CH4_t, ECH4)
        plt.savefig('results/plot-CH4.svg')
    if plot == 'all' or plot == 'N2O' or 'N2O' in plot:
        plot_N2O(D_N2O_t, D_HVSNK_N2O_t, D_EBB_N2O_t, EN2O)
        plt.savefig('results/plot-N20.svg')
    if plot == 'all' or plot == 'O3' or 'O3' in plot:
        plot_O3(D_O3t_t, D_O3s_t, D_EESC_t, D_N2O_lag_t, D_gst_t)
        plt.savefig('results/plot-03.svg')
    if plot == 'all' or plot == 'AER' or 'AER' in plot:
        plot_AER(D_SO4_t, D_POA_t, D_BC_t, D_NO3_t, D_SOA_t, D_AERh_t, RF_SO4_t, RF_POA_t,
                 RF_BC_t, RF_NO3_t, RF_SOA_t, RF_cloud_t)
        plt.savefig('results/plot-AER.svg')
    if plot == 'all' or plot == 'clim' or 'clim' in plot:
        plot_clim(RF_t, D_gst_t, D_gyp_t, RF_CO2_t, RF_CH4_t, RF_H2Os_t, RF_N2O_t,
                  RF_halo_t, RF_O3t_t, RF_O3s_t, RF_SO4_t, RF_POA_t, RF_BC_t, RF_NO3_t,
                  RF_SOA_t, RF_cloud_t, RF_BCsnow_t, RF_LCC_t, RFcon, RFvolc, RFsolar)
        plt.savefig('results/plot-clim.svg')

    # ===========
    # D. OUTPUTS
    # ===========

    output = []
    for var in var_output:
        exec('output.append(' + var + '_t)')

    return output


##################################################
#   2. CONTROL PLOTS
##################################################

# =========
# 2.1. CO2
# =========

def plot_CO2(D_CO2, OSNK, LSNK, ELUC, EFF, D_AREA, D_npp, D_efire, D_fmort, D_rh1, D_fmet,
             D_rh2, D_FIN, D_FOUT, D_FCIRC, D_MORT_luc, D_EFIRE_luc, D_RH1_luc, D_RH2_luc,
             EHWP1_luc, EHWP2_luc, EHWP3_luc):
    plt.figure()

    # atmospheric CO2
    ax = plt.subplot(2, 3, 1)
    plt.plot(1700 + np.arange(ind_final + 1), D_CO2, color='k', lw=2, label='OSCAR')
    plt.plot(1700 + np.arange(len(CO2_ipcc)), CO2_ipcc - CO2_0, color='r', lw=2, ls='--',
             label='IPCC')
    if (ind_final > ind_cdiac):
        plt.plot(1700 + np.arange(min(len(CO2_rcp), ind_final + 1)),
                 CO2_rcp[:min(len(CO2_rcp), ind_final + 1), 0] - CO2_0, color='0.8', lw=2,
                 ls=':', label='RCP2.6')
        plt.plot(1700 + np.arange(min(len(CO2_rcp), ind_final + 1)),
                 CO2_rcp[:min(len(CO2_rcp), ind_final + 1), 1] - CO2_0, color='0.6', lw=2,
                 ls=':', label='RCP4.5')
        plt.plot(1700 + np.arange(min(len(CO2_rcp), ind_final + 1)),
                 CO2_rcp[:min(len(CO2_rcp), ind_final + 1), 2] - CO2_0, color='0.4', lw=2,
                 ls=':', label='RCP6.0')
        plt.plot(1700 + np.arange(min(len(CO2_rcp), ind_final + 1)),
                 CO2_rcp[:min(len(CO2_rcp), ind_final + 1), 3] - CO2_0, color='0.2', lw=2,
                 ls=':', label='RCP8.5')
    plt.title('$\Delta$CO2 (ppm)', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + ind_final])

    # budget fluxes
    ax = plt.subplot(2, 3, 2)
    plt.plot([1700, 1700 + ind_final + 1], [0, 0], 'k-')
    plt.plot(1700 + np.arange(ind_final + 1), np.sum(EFF, 1), color='#666666', lw=2,
             label='EFF')
    plt.plot(1700 + np.arange(ind_final + 1), np.sum(ELUC, 1), color='#993300', lw=2,
             label='ELUC')
    plt.plot(1700 + np.arange(ind_final + 1), OSNK, color='#000099', lw=2, label='OSNK')
    plt.plot(1700 + np.arange(ind_final + 1), LSNK, color='#009900', lw=2, label='LSNK')
    plt.plot(1700 + np.arange(ind_final) + 1, alpha_CO2 * (D_CO2[1:] - D_CO2[:-1]),
             color='#FFCC00', lw=2, label='d_CO2')
    plt.plot(1700 + np.arange(len(EFF_gcp)), EFF_gcp, color='#666666', ls='--')
    plt.plot(1700 + np.arange(len(ELUC_gcp)), ELUC_gcp, color='#CC3300', ls='--')
    plt.plot(1700 + np.arange(len(OSNK_gcp)), OSNK_gcp, color='#000099', ls='--')
    plt.plot(1700 + np.arange(len(LSNK_gcp)), LSNK_gcp, color='#009900', ls='--')
    plt.plot(1700 + np.arange(len(d_CO2_gcp)), d_CO2_gcp, color='#FFCC00', ls='--')
    plt.plot([1700, 1700], [0, 0], 'k--', label='GCP')
    plt.title('CO2 fluxes (GtC/yr)', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + ind_final])

    # airborne fraction
    ax = plt.subplot(2, 3, 3)
    plt.plot([1700, 1700 + ind_final + 1], [0, 0], 'k-')
    plt.plot(1700 + np.arange(ind_final) + 1,
             alpha_CO2 * (D_CO2[1:] - D_CO2[:-1]) / np.sum(EFF + ELUC, 1)[1:],
             color='#FFCC00', lw=1, label='AF')
    plt.plot(1700 + np.arange(ind_final + 1), -OSNK / np.sum(EFF + ELUC, 1),
             color='#000099', lw=1, label='OF')
    plt.plot(1700 + np.arange(ind_final + 1), -LSNK / np.sum(EFF + ELUC, 1),
             color='#009900', lw=1, label='LF')
    plt.plot(np.arange(1959, 1700 + ind_cdiac + 1),
             np.ones([ind_cdiac - 259 + 1]) * np.mean(
                 (alpha_CO2 * (D_CO2[1:] - D_CO2[:-1]) / np.sum(EFF + ELUC, 1)[1:])[
                 259 - 1:ind_cdiac]), color='k', lw=2, label='OSCAR')
    plt.plot(np.arange(1959, 1700 + ind_cdiac + 1),
             np.ones([ind_cdiac - 259 + 1]) * np.mean(
                 (d_CO2_gcp / (EFF_gcp + ELUC_gcp))[259:ind_cdiac + 1]), color='r', lw=2,
             ls='--', label='GCP')
    plt.title('airborne fraction', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + ind_final])
    ax.set_ylim([-0.2, 1.2])

    # ELUC details
    ax = plt.subplot(2, 3, 4)
    plt.plot([1700, 1700 + ind_final + 1], [0, 0], 'k-')
    plt.plot(1700 + np.arange(ind_final + 1), np.sum(ELUC, 1), color='k', ls='-.', lw=2,
             label='ELUC')
    plt.plot(1700 + np.arange(ind_final + 1),
             np.sum(np.sum(np.sum(np.sum(D_EFIRE_luc + D_RH1_luc + D_RH2_luc, 4), 3), 2),
                    1), color='#009900', lw=2, label='ELUC_bio')
    plt.plot(1700 + np.arange(ind_final + 1),
             np.sum(np.sum(np.sum(np.sum(EHWP1_luc + EHWP2_luc + EHWP3_luc, 4), 3), 2),
                    1), color='#993300', lw=2, label='ELUC_hwp')
    plt.plot(1700 + np.arange(ind_final + 1),
             np.sum(np.sum(np.sum(np.sum(EHWP1_luc, 4), 3), 2), 1), color='#FF3300', lw=1,
             label='EHWP1')
    plt.plot(1700 + np.arange(ind_final + 1),
             np.sum(np.sum(np.sum(np.sum(EHWP2_luc, 4), 3), 2), 1), color='#CC9900', lw=1,
             label='EHWP2')
    plt.plot(1700 + np.arange(ind_final + 1),
             np.sum(np.sum(np.sum(np.sum(EHWP3_luc, 4), 3), 2), 1), color='#663300', lw=1,
             label='EHWP3')
    plt.title('ELUC fluxes (GtC/yr)', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))

    # LSNK details
    ax = plt.subplot(2, 3, 5)
    plt.plot([1700, 1700 + ind_final + 1], [0, 0], 'k-')
    plt.plot(1700 + np.arange(ind_final + 1), -LSNK, color='k', lw=2, ls='-.',
             label='$-$LSNK')
    plt.plot(1700 + np.arange(ind_final + 1),
             np.sum(np.sum(D_npp * (AREA_0 + D_AREA), 2), 1), color='#009900', lw=2,
             label='D_NPP')
    plt.plot(1700 + np.arange(ind_final + 1),
             np.sum(np.sum(D_efire * (AREA_0 + D_AREA), 2), 1), color='#FF3300', lw=2,
             label='D_EFIRE')
    plt.plot(1700 + np.arange(ind_final + 1),
             np.sum(np.sum(D_fmort * (AREA_0 + D_AREA), 2), 1), color='#336633', lw=2,
             label='D_FMORT')
    plt.plot(1700 + np.arange(ind_final + 1),
             np.sum(np.sum((D_rh1 + D_rh2) * (AREA_0 + D_AREA), 2), 1), color='#663300',
             lw=2, label='D_RH')
    plt.title('LSNK fluxes (GtC/yr)', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))

    # OSNK details
    ax = plt.subplot(2, 3, 6)
    plt.plot([1700, 1700 + ind_final + 1], [0, 0], 'k-')
    plt.plot(1700 + np.arange(ind_final + 1), -OSNK, color='k', lw=2, ls='-.',
             label='$-$OSNK')
    plt.plot(1700 + np.arange(ind_final + 1), np.sum(D_FIN, 1), color='#000099', lw=2,
             label='D_FIN')
    plt.plot(1700 + np.arange(ind_final + 1), np.sum(D_FOUT, 1), color='#0099FF', lw=2,
             label='D_FOUT')
    plt.plot(1700 + np.arange(ind_final + 1), np.sum(D_FCIRC, 1), color='#663399', lw=2,
             label='D_FCIRC')
    plt.title('OSNK fluxes (GtC/yr)', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))


# =========
# 2.2. CH4
# =========

def plot_CH4(D_CH4, D_OHSNK_CH4, D_HVSNK_CH4, D_XSNK_CH4, D_EWET, D_EBB_CH4, ECH4):
    plt.figure()

    # atmospheric CH4
    ax = plt.subplot(2, 3, 1)
    plt.plot(1700 + np.arange(ind_final + 1), D_CH4, color='k', lw=2, label='OSCAR')
    plt.plot(1700 + np.arange(len(CH4_ipcc)), CH4_ipcc - CH4_0, color='r', lw=2, ls='--',
             label='IPCC')
    if (ind_final > ind_cdiac):
        plt.plot(1700 + np.arange(min(len(CH4_rcp), ind_final + 1)),
                 CH4_rcp[:min(len(CH4_rcp), ind_final + 1), 0] - CH4_0, color='0.8', lw=2,
                 ls=':', label='RCP2.6')
        plt.plot(1700 + np.arange(min(len(CH4_rcp), ind_final + 1)),
                 CH4_rcp[:min(len(CH4_rcp), ind_final + 1), 1] - CH4_0, color='0.6', lw=2,
                 ls=':', label='RCP4.5')
        plt.plot(1700 + np.arange(min(len(CH4_rcp), ind_final + 1)),
                 CH4_rcp[:min(len(CH4_rcp), ind_final + 1), 2] - CH4_0, color='0.4', lw=2,
                 ls=':', label='RCP6.0')
        plt.plot(1700 + np.arange(min(len(CH4_rcp), ind_final + 1)),
                 CH4_rcp[:min(len(CH4_rcp), ind_final + 1), 3] - CH4_0, color='0.2', lw=2,
                 ls=':', label='RCP8.5')
    plt.title('$\Delta$CH4 (ppb)', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + ind_final])

    # budget fluxes
    ax = plt.subplot(2, 3, 2)
    plt.plot([1700, 1700 + ind_final + 1], [0, 0], 'k-')
    plt.plot(1700 + np.arange(ind_final + 1), np.sum(ECH4, 1), color='#666666', lw=2,
             label='ECH4')
    plt.plot(1700 + np.arange(ind_final + 1), np.sum(D_EBB_CH4, 1), color='#993300', lw=2,
             label='D_EBB')
    plt.plot(1700 + np.arange(ind_final + 1), np.sum(D_EWET, 1), color='#006666', lw=2,
             label='D_EWET')
    plt.plot(1700 + np.arange(ind_final + 1), (D_OHSNK_CH4 + D_HVSNK_CH4 + D_XSNK_CH4),
             color='#990066', lw=2, label='D_SNK')
    plt.plot(1700 + np.arange(ind_final) + 1, alpha_CH4 * (D_CH4[1:] - D_CH4[:-1]),
             color='#FFCC00', lw=2, label='d_CH4')
    plt.title('CH4 fluxes (MtC/yr)', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + ind_final])

    # lifetime
    ax = plt.subplot(2, 3, 3)
    plt.plot(1700 + np.arange(ind_final + 1), alpha_CH4 * (CH4_0 + D_CH4) / (
                alpha_CH4 * CH4_0 * (
                    1 / tau_CH4_OH + 1 / tau_CH4_hv + 1 / tau_CH4_soil + 1 / tau_CH4_ocean) - D_OHSNK_CH4 - D_HVSNK_CH4 - D_XSNK_CH4),
             color='k', lw=2, label='OSCAR')
    plt.plot(1700 + np.arange(ind_final + 1),
             alpha_CH4 * (CH4_0 + D_CH4) / (alpha_CH4 * CH4_0 / tau_CH4_OH - D_OHSNK_CH4),
             color='k', lw=1, label='OH only')
    plt.title('CH4 lifetime (yr)', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + ind_final])

    # wetlands
    ax = plt.subplot(2, 3, 4)
    plt.title('wetlands', fontsize='medium')

    # biomass burning
    ax = plt.subplot(2, 3, 5)
    plt.title('biomass burning', fontsize='medium')


# =========
# 2.3. N2O
# =========

def plot_N2O(D_N2O, D_HVSNK_N2O, D_EBB_N2O, EN2O):
    plt.figure()

    # atmospheric N2O
    ax = plt.subplot(2, 3, 1)
    plt.plot(1700 + np.arange(ind_final + 1), D_N2O, color='k', lw=2, label='OSCAR')
    plt.plot(1700 + np.arange(len(N2O_ipcc)), N2O_ipcc - N2O_0, color='r', lw=2, ls='--',
             label='IPCC')
    if (ind_final > ind_cdiac):
        plt.plot(1700 + np.arange(min(len(N2O_rcp), ind_final + 1)),
                 N2O_rcp[:min(len(N2O_rcp), ind_final + 1), 0] - N2O_0, color='0.8', lw=2,
                 ls=':', label='RCP2.6')
        plt.plot(1700 + np.arange(min(len(N2O_rcp), ind_final + 1)),
                 N2O_rcp[:min(len(N2O_rcp), ind_final + 1), 1] - N2O_0, color='0.6', lw=2,
                 ls=':', label='RCP4.5')
        plt.plot(1700 + np.arange(min(len(N2O_rcp), ind_final + 1)),
                 N2O_rcp[:min(len(N2O_rcp), ind_final + 1), 2] - N2O_0, color='0.4', lw=2,
                 ls=':', label='RCP6.0')
        plt.plot(1700 + np.arange(min(len(N2O_rcp), ind_final + 1)),
                 N2O_rcp[:min(len(N2O_rcp), ind_final + 1), 3] - N2O_0, color='0.2', lw=2,
                 ls=':', label='RCP8.5')
    plt.title('$\Delta$N2O (ppb)', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + ind_final])

    # budget fluxes
    ax = plt.subplot(2, 3, 2)
    plt.plot([1700, 1700 + ind_final + 1], [0, 0], 'k-')
    plt.plot(1700 + np.arange(ind_final + 1), np.sum(EN2O, 1), color='#666666', lw=2,
             label='EN2O')
    plt.plot(1700 + np.arange(ind_final + 1), np.sum(D_EBB_N2O, 1), color='#993300', lw=2,
             label='D_EBB')
    plt.plot(1700 + np.arange(ind_final + 1), D_HVSNK_N2O, color='#990066', lw=2,
             label='D_SNK')
    plt.plot(1700 + np.arange(ind_final) + 1, alpha_N2O * (D_N2O[1:] - D_N2O[:-1]),
             color='#FFCC00', lw=2, label='d_N2O')
    plt.title('N2O fluxes (MtN/yr)', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + ind_final])

    # lifetime
    ax = plt.subplot(2, 3, 3)
    plt.plot(1700 + np.arange(ind_final + 1),
             alpha_N2O * (N2O_0 + D_N2O) / (alpha_N2O * N2O_0 / tau_N2O_hv - D_HVSNK_N2O),
             color='k', lw=2, label='OSCAR')
    plt.title('N2O lifetime (yr)', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + ind_final])


# ========
# 2.4. O3
# ========

def plot_O3(D_O3t, D_O3s, D_EESC, D_N2O_lag, D_gst):
    plt.figure()

    # tropospheric O3
    ax = plt.subplot(2, 3, 1)
    plt.plot(1700 + np.arange(ind_final + 1), D_O3t, color='k', lw=2, label='OSCAR')
    plt.title('$\Delta$O3 trop. (DU)', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + ind_final])

    # stratospheric O3
    ax = plt.subplot(2, 3, 2)
    plt.plot(1700 + np.arange(ind_final + 1), D_O3s, color='k', lw=2, label='OSCAR')
    plt.title('$\Delta$O3 strat. (DU)', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + ind_final])

    # EESC
    ax = plt.subplot(2, 3, 3)
    plt.plot(1700 + np.arange(ind_final + 1), D_EESC, color='k', lw=2, label='OSCAR')
    plt.plot(1700 + np.arange(ind_final + 1),
             (chi_O3s_N2O * D_N2O_lag * (1 - D_EESC / EESC_x) / chi_O3s_EESC), color='k',
             lw=1, label='N2O effect')
    plt.title('$\Delta$EESC (ppt)', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + ind_final])

    # age-of-air
    ax = plt.subplot(2, 3, 4)
    plt.plot(1700 + np.arange(ind_final + 1), tau_lag / (1 + gamma_age * D_gst),
             color='k', lw=2, label='OSCAR')
    plt.title('mean age-of-air (yr)', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + ind_final])


# ==============
# 2.5. Aerosols
# ==============

def plot_AER(D_SO4, D_POA, D_BC, D_NO3, D_SOA, D_AERh, RF_SO4, RF_POA, RF_BC, RF_NO3,
             RF_SOA, RF_cloud):
    plt.figure()

    # atmospheric burden
    ax = plt.subplot(2, 3, 1)
    plt.plot(1700 + np.arange(ind_final + 1), D_SO4, color='b', lw=2, label='D_SO4')
    plt.plot(1700 + np.arange(ind_final + 1), D_POA, color='m', lw=2, label='D_POA')
    plt.plot(1700 + np.arange(ind_final + 1), D_BC, color='r', lw=2, label='D_BC')
    plt.plot(1700 + np.arange(ind_final + 1), D_NO3, color='g', lw=2, label='D_NO3')
    plt.plot(1700 + np.arange(ind_final + 1), D_SOA, color='y', lw=2, label='D_SOA')
    plt.plot(1700 + np.arange(ind_final + 1), D_AERh, color='c', lw=2, label='D_AERh')
    plt.title('$\Delta$ burdens (Tg)', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + ind_final])

    # radiative forcing
    ax = plt.subplot(2, 3, 4)
    plt.plot([1700, 1700 + ind_final], [0, 0], 'k-')
    plt.plot(1700 + np.arange(ind_final + 1), RF_SO4, color='b', lw=2, label='RF_SO4')
    plt.errorbar([2010], [-0.40], yerr=[[0.20], [0.20]], marker='o', mfc='b', color='k')
    plt.plot(1700 + np.arange(ind_final + 1), RF_POA, color='m', lw=2, label='RF_POA')
    plt.errorbar([2010], [-0.29], yerr=[[-0.29 * 0.63], [-0.29 * 0.72]], marker='o',
                 mfc='m', color='k')
    plt.plot(1700 + np.arange(ind_final + 1), RF_BC, color='r', lw=2, label='RF_BC')
    plt.errorbar([2010], [+0.60], yerr=[[+0.60 * 0.61], [+0.60 * 0.70]], marker='o',
                 mfc='r', color='k')
    plt.plot(1700 + np.arange(ind_final + 1), RF_NO3, color='g', lw=2, label='RF_NO3')
    plt.errorbar([2010], [-0.11], yerr=[[0.19], [0.08]], marker='o', mfc='g', color='k')
    plt.plot(1700 + np.arange(ind_final + 1), RF_SOA, color='y', lw=2, label='RF_SOA')
    plt.errorbar([2010], [-0.03], yerr=[[0.24], [0.23]], marker='o', mfc='y', color='k')
    plt.plot(1700 + np.arange(ind_final + 1), RF_cloud, color='c', lw=2, label='RF_cloud')
    plt.errorbar([2010], [-0.45], yerr=[[0.75], [0.45]], marker='o', mfc='c', color='k')
    # plt.errorbar([2010],[-0.10],yerr=[[0.20],[0.20]],marker='o',mfc='0.5',color='k')
    plt.title('RF (W/m2)', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, max(1700 + ind_final, 2010 + 10)])


# =============
# 2.6. Climate
# =============

def plot_clim(RF, D_gst, D_gyp, RF_CO2, RF_CH4, RF_H2Os, RF_N2O, RF_halo, RF_O3t, RF_O3s,
              RF_SO4, RF_POA, RF_BC, RF_NO3, RF_SOA, RF_cloud, RF_BCsnow, RF_LCC, RFcon,
              RFvolc, RFsolar):
    plt.figure()

    # radiative forcing
    ax = plt.subplot(2, 3, 1)
    plt.plot([1700, 1700 + ind_final], [0, 0], 'k-')
    plt.plot(1700 + np.arange(ind_final + 1), RF, color='k', lw=2, label='OSCAR')
    plt.plot(1700 + np.arange(len(RF_ipcc)), RF_ipcc, color='r', lw=2, ls='--',
             label='IPCC')
    plt.title('RF (W/m2)', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + ind_final])

    # global temperature
    ax = plt.subplot(2, 3, 2)
    plt.plot([1700, 1700 + ind_final], [0, 0], 'k-')
    plt.plot(1700 + np.arange(ind_final + 1), D_gst - np.mean(D_gst[200:230]), color='k',
             lw=2, label='OSCAR')
    plt.plot(1700 + np.arange(len(gst_giss)), gst_giss - np.mean(gst_giss[200:230]),
             color='b', ls='--', label='GISS')
    plt.plot(1700 + np.arange(len(gst_had)), gst_had - np.mean(gst_had[200:230]),
             color='g', ls='--', label='Hadley')
    plt.plot(1700 + np.arange(len(gst_ncdc)), gst_ncdc - np.mean(gst_ncdc[200:230]),
             color='m', ls='--', label='NCDC')
    plt.title('$\Delta$ temp. (K)', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + ind_final])

    # global precipitations
    ax = plt.subplot(2, 3, 3)
    plt.plot(1700 + np.arange(ind_final + 1), D_gyp, color='k', lw=2, label='OSCAR')
    plt.title('$\Delta$ precip. (mm)', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + ind_final])

    # RF details
    ax = plt.subplot(2, 3, 4)
    plt.plot([1700, 1700 + ind_final], [0, 0], 'k-')
    plt.plot(1700 + np.arange(ind_final + 1),
             RF_CO2 + RF_CH4 + RF_N2O + RF_halo + RF_H2Os, color='r', lw=2, label='WMGHG')
    plt.plot(1700 + np.arange(len(RF_WMGHG_ipcc)), RF_WMGHG_ipcc, color='r', ls='--')
    plt.plot(1700 + np.arange(ind_final + 1), RF_O3t + RF_O3s, color='y', lw=2,
             label='O3')
    plt.plot(1700 + np.arange(len(RF_O3_ipcc)), RF_O3_ipcc, color='y', ls='--')
    plt.plot(1700 + np.arange(ind_final + 1),
             RF_SO4 + RF_POA + RF_BC + RF_NO3 + RF_SOA + RF_cloud, color='b', lw=2,
             label='AER')
    plt.plot(1700 + np.arange(len(RF_AER_ipcc)), RF_AER_ipcc, color='b', ls='--')
    plt.plot(1700 + np.arange(ind_final + 1), RF_BCsnow + RF_LCC, color='g', lw=2,
             label='Alb.')
    plt.plot(1700 + np.arange(len(RF_Alb_ipcc)), RF_Alb_ipcc, color='g', ls='--')
    plt.plot(1700 + np.arange(ind_final + 1), RFcon, color='k', ls='--', label='Ant.')
    plt.plot(1700 + np.arange(ind_final + 1), RFvolc + RFsolar, color='0.5', ls='--',
             label='Nat.')
    plt.title('RF (W/m2)', fontsize='medium')
    plt.legend(loc=0, ncol=2, prop=FontProperties(size='small'))
    plt.xticks(rotation=27)
    ax.set_xlim([1700, 1700 + ind_final])

