'''
Created on Oct 31, 2021

@author: mzwier
'''

import pandas as pd

_elements_table = [
    (1,'H','Hydrogen',1.008,1),
    (2,'He','Helium',4.0026,4),
    (3,'Li','Lithium',6.94,7),
    (4,'Be','Beryllium',9.0122,9),
    (5,'B','Boron',10.81,11),
    (6,'C','Carbon',12.011,12),
    (7,'N','Nitrogen',14.007,14),
    (8,'O','Oxygen',15.999,16),
    (9,'F','Fluorine',18.998,19),
    (10,'Ne','Neon',20.18,20),
    (11,'Na','Sodium',22.99,23),
    (12,'Mg','Magnesium',24.305,24),
    (13,'Al','Aluminium',26.982,27),
    (14,'Si','Silicon',28.085,28),
    (15,'P','Phosphorus',30.974,31),
    (16,'S','Sulfur',32.06,32),
    (17,'Cl','Chlorine',35.45,35),
    (18,'Ar','Argon',39.95,40),
    (19,'K','Potassium',39.098,93),
    (20,'Ca','Calcium',40.078,40),
    (21,'Sc','Scandium',44.956,45),
    (22,'Ti','Titanium',47.867,48),
    (23,'V','Vanadium',50.942,51),
    (24,'Cr','Chromium',51.996,52),
    (25,'Mn','Manganese',54.938,55),
    (26,'Fe','Iron',55.845,56),
    (27,'Co','Cobalt',58.933,59),
    (28,'Ni','Nickel',58.693,58),
    (29,'Cu','Copper',63.546,63),
    (30,'Zn','Zinc',65.38,64),
    (31,'Ga','Gallium',69.723,69),
    (32,'Ge','Germanium',72.63,74),
    (33,'As','Arsenic',74.922,75),
    (34,'Se','Selenium',78.971,80),
    (35,'Br','Bromine',79.904,79),
    (36,'Kr','Krypton',83.798,84),
    (37,'Rb','Rubidium',85.468,85),
    (38,'Sr','Strontium',87.62,88),
    (39,'Y','Yttrium',88.906,89),
    (40,'Zr','Zirconium',91.224,90),
    (41,'Nb','Niobium',92.906,93),
    (42,'Mo','Molybdenum',95.95,98),
    (43,'Tc','Technetium',97.0,97),
    (44,'Ru','Ruthenium',101.07,102),
    (45,'Rh','Rhodium',102.91,103),
    (46,'Pd','Palladium',106.42,106),
    (47,'Ag','Silver',107.87,107),
    (48,'Cd','Cadmium',112.41,114),
    (49,'In','Indium',114.82,115),
    (50,'Sn','Tin',118.71,120),
    (51,'Sb','Antimony',121.76,121),
    (52,'Te','Tellurium',127.6,130),
    (53,'I','Iodine',126.9,127),
    (54,'Xe','Xenon',131.29,132),
    (55,'Cs','Caesium',132.91,133),
    (56,'Ba','Barium',137.33,138),
    (57,'La','Lanthanum',138.91,139),
    (58,'Ce','Cerium',140.12,140),
    (59,'Pr','Praseodymium',140.91,141),
    (60,'Nd','Neodymium',144.24,142),
    (61,'Pm','Promethium',145.0,145),
    (62,'Sm','Samarium',150.36,152),
    (63,'Eu','Europium',151.96,153),
    (64,'Gd','Gadolinium',157.25,158),
    (65,'Tb','Terbium',158.93,159),
    (66,'Dy','Dysprosium',162.5,164),
    (67,'Ho','Holmium',164.93,165),
    (68,'Er','Erbium',167.26,166),
    (69,'Tm','Thulium',168.93,169),
    (70,'Yb','Ytterbium',173.05,174),
    (71,'Lu','Lutetium',174.97,175),
    (72,'Hf','Hafnium',178.49,180),
    (73,'Ta','Tantalum',180.95,181),
    (74,'W','Tungsten',183.84,184),
    (75,'Re','Rhenium',186.21,187),
    (76,'Os','Osmium',190.23,192),
    (77,'Ir','Iridium',192.22,193),
    (78,'Pt','Platinum',195.08,195),
    (79,'Au','Gold',196.97,197),
    (80,'Hg','Mercury',200.59,202),
    (81,'Tl','Thallium',204.38,205),
    (82,'Pb','Lead',207.2,206),
    (83,'Bi','Bismuth',208.98,209),
    (84,'Po','Polonium',209.0,209),
    (85,'At','Astatine',210.0,210),
    (86,'Rn','Radon',222.0,222),
    (87,'Fr','Francium',223.0,223),
    (88,'Ra','Radium',226.0,226),
    (89,'Ac','Actinium',227.0,227),
    (90,'Th','Thorium',232.04,232),
    (91,'Pa','Protactinium',231.04,231),
    (92,'U','Uranium',238.03,238),
    (93,'Np','Neptunium',237.0,237),
    (94,'Pu','Plutonium',244.0,244),
    (95,'Am','Americium',243.0,243),
    (96,'Cm','Curium',247.0,247),
    (97,'Bk','Berkelium',247.0,247),
    (98,'Cf','Californium',251.0,251),
    (99,'Es','Einsteinium',252.0,252),
    (100,'Fm','Fermium',257.0,257),
    (101,'Md','Mendelevium',258.0,258),
    (102,'No','Nobelium',259.0,259),
    (103,'Lr','Lawrencium',266.0,266),
    (104,'Rf','Rutherfordium',267.0,267),
    (105,'Db','Dubnium',268.0,268),
    (106,'Sg','Seaborgium',269.0,269),
    (107,'Bh','Bohrium',270.0,270),
    (108,'Hs','Hassium',269.0,269),
    (109,'Mt','Meitnerium',278.0,278),
    (110,'Ds','Darmstadtium',281.0,281),
    (111,'Rg','Roentgenium',282.0,282),
    (112,'Cn','Copernicium',285.0,285),
    (113,'Nh','Nihonium',286.0,286),
    (114,'Fl','Flerovium',289.0,289),
    (115,'Mc','Moscovium',290.0,290),
    (116,'Lv','Livermorium',293.0,293),
    (117,'Ts','Tennessine',294.0,294),
    (118,'Og','Oganesson',294.0,294),
]

elements = pd.DataFrame.from_records(_elements_table, 
                                     index='number', 
                                     columns=['number', 'symbol', 'name', 'mass', 'prevalent_isotope' ])

elements_by_number = elements
elements_by_symbol = elements.reset_index().set_index('symbol')
elements_by_name   = elements.reset_index().set_index('name')

def element_symbol_to_number(sym):
    sym = sym.title()
    return elements_by_symbol['number'][sym]

def element_number_to_symbol(num):
    num = int(num)
    return elements_by_number['symbol'][num]

def element_label_to_number(label):
    '''Attempt to convert a label (name, symbol, or atomic number) to atomic number'''
    
    # First, try some fast heuristics
    # If any of them fail, be more careful
    
    try:
        if isinstance(label, str):
            if len(label) <= 3:
                # Probably a symbol
                return elements_by_symbol['number'][label.title()]
            else:
                # Probaably a name
                return elements_by_name['number'][label.title()]
        else:
            # Probably an atomic number
            return int(label)
    except (IndexError,ValueError):
        pass
    
    # try atomic symbol
    try:
        return elements_by_symbol['number'][label.title()]
    except KeyError:
        pass
    
    # Try atomic number
    try:
        return int(label)
    except ValueError:
        pass
    
    # Try name
    try:
        return elements_by_name['number'][label.title()]
    except KeyError:
        pass
    
    raise ValueError('could not identify {!r} as an element'.format(label))
        
        
    
