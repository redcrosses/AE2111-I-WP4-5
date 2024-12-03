from WP4_1.main1 import main1
from WP4_2.main2 import main2
from WP4_3.main3 import main3

n_positive, n_negative = main3()
# print(n_positive)
loads_positive, loads_negative, spanwise_position = main1(n_positive, n_negative)
main2(loads_positive, spanwise_position)




