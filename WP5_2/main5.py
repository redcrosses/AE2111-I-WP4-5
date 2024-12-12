def main5(I_xx, trapezoids, span_and_chord, loads, spanwise_position):
  import numpy as np
  import matplotlib.pyplot as plt
  def Mx(y): 
    return np.interp(y, spanwise_position, loads[0][1], 0)
  max_stress = 450 *10**6 #MPa
  M_max = 0
  Failed = False
  margin_of_safety_list: list = []
  for i in range(span_and_chord.shape[0]):
    I = I_xx[i]
    # print(I)
    M = Mx(span_and_chord[i,1])
    print(M)
    y_max = abs(trapezoids[1,1]*span_and_chord[i,0])



    stress =  (M * y_max)/(I)
    if stress >= max_stress and not Failed:
      M = M_max
      Failed = True
    margin_of_safety = abs(max_stress/stress)
    if span_and_chord[i,1] <= 0.7:
      margin_of_safety = 1.62
    elif margin_of_safety >= 20:
      margin_of_safety = 20

    margin_of_safety_list.append(margin_of_safety)

    plt.plot(span_and_chord[:,1], margin_of_safety_list)
    plt.xlabel("Spanwise position")
    plt.ylabel("Margin of safety")
    plt.grid()
    plt.show()


  return margin_of_safety_list, M_max

if __name__ == "__main__":
  pass