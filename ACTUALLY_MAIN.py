from WP4_1.main1 import main1
from WP4_2.main2 import main2
from WP4_3.main3 import main3
import sys

if not sys.warnoptions:
    import warnings
    warnings.simplefilter("ignore")

print("MAIN3: Finding max, min loads...")
n_positive, n_negative = main3()
# print(n_positive)

print("MAIN1: Finding planform loading diagrams...")
loads_positive, loads_negative, spanwise_position = main1(n_positive, n_negative)
# print(loads_positive, loads_negative)
print("MAIN2:")
main2((loads_positive, loads_negative), spanwise_position, (n_positive, n_negative))
