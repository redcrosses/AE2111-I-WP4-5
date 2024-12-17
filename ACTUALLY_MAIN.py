from WP4_1.main1 import main1
from WP4_2.main2 import main2
from WP4_3.main3 import main3
from WP5_1.main4 import main4
from WP5_2.main5 import main5
import sys

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

print("\033[96m {} \033[00m".format("WP4.3:"))
print("Finding load factors")
n_positive, n_negative = main3()
print(n_positive)
print(n_negative)

print("Finding planform loading diagrams...")
loads_positive, loads_negative, spanwise_position = main1(n_positive, n_negative)
# print(loads_positive, loads_negative)
design: object = main2((loads_positive, loads_negative), spanwise_position, (n_positive, n_negative), 0.1056927, 0.07702, 0.005, 0.01,  42, 2e-4)

margin, max = main5(design.moi_x_list, design.trapezoid,  design.chords_along_span, (loads_positive, loads_negative),spanwise_position)
main4(design.moi_x_list, design.trapezoid, design.stringers, design.chords_along_span, (loads_positive, loads_negative), spanwise_position, design)

# moi_x_list, trapezoid, stringer_positions, span_positions_and_chord 
# moi_x_list, trapezoid, stringer_positions, span_positions_and_chord = main2((loads_positive, loads_negative), spanwise_position, (n_positive, n_negative), 0.12079, 0.07702, 0.01, 0.01, 42, 4e-4)
# margin, max = main5(moi_x_list, trapezoid, span_positions_and_chord, (loads_positive, loads_negative),spanwise_position)

# moi_x_list, trapezoid, stringer_positions, span_positions_and_chord = main2((loads_positive, loads_negative), spanwise_position, (n_positive, n_negative), 0.12079, 0.06808, 0.01, 0.015, 42, 4e-4)
# margin, max = main5(moi_x_list, trapezoid, span_positions_and_chord, (loads_positive, loads_negative),spanwise_position)
# main2((loads_positive, loads_negative), spanwise_position, (n_positive, n_negative), 0.12079, 0.07702, 0.01, 0.01, 42, 4e-4)
# main2((loads_positive, loads_negative), spanwise_position, (n_positive, n_negative), 0.12079, 0.06808, 0.01, 0.015, 42, 4e-4)
# margin, max = main5(moi_x_list, trapezoid, span_positions_and_chord, (loads_positive, loads_negative),spanwise_position)
# print(margin)