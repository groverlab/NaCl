import Tkinter
import numpy
import scipy.optimize

root = Tkinter.Tk()

choice = Tkinter.StringVar(root)
choice.set("select...")

def callback():
    print "clicked", solute_mass_entry.get()

def o_changed(throwaway):
    print "calculating for", choice.get(), "..."
    if choice.get() == "NaCl":
        temperatures = [20.0, 25.0, 30.0, 40.0]
        molalities = numpy.array([0.100, 0.250, 0.500, 0.750, 1.000, 2.000, 3.000, 4.000, 5.000])
        specific_volumes = numpy.array([#[0.995732, 0.989259, 0.978889, 0.968991, 0.959525, 0.925426, 0.896292, 0.870996, 0.848646],  # 0.0 C
                                        #[0.995998, 0.989781, 0.979804, 0.970256, 0.961101, 0.927905, 0.899262, 0.874201, 0.851958],  # 10.0 C
                                        [0.997620, 0.991564, 0.981833, 0.972505, 0.963544, 0.930909, 0.902565, 0.877643, 0.855469],  # 20.0 C
                                        [0.998834, 0.992832, 0.983185, 0.973932, 0.965038, 0.932590, 0.904339, 0.879457, 0.857301],  # 25.0 C
                                        [1.000279, 0.994319, 0.984735, 0.975539, 0.966694, 0.934382, 0.906194, 0.881334, 0.859185],  # 30.0 C
                                        [1.003796, 0.997883, 0.988374, 0.979243, 0.970455, 0.938287, 0.910145, 0.885276, 0.863108]]) # 40.0 C

        densities = 1.0/specific_volumes 
        mass_NaCl_over_mass_water = molalities * 58.443 / 1000.0

        def error_function((a, b, c, d, e)):
            sum = 0.0
            for t, i in zip(temperatures, range(len(temperatures))):
                for m, j in zip(mass_NaCl_over_mass_water, range(len(mass_NaCl_over_mass_water))):
                    sum += ((a*t**2 + b*t + c*m**2 + d*m + e) - densities[i][j])**2
            return sum

        a, b, c, d, e = scipy.optimize.fmin(error_function, (0.0001, 2, 3, 4, 1), full_output=1, xtol=1e-29, ftol=1e-29, maxiter=10000, maxfun=10000)[0]

        x = []
        y = []
        ztable = []
        zfunc = []
        count = 0
        max_diff = 0.0
        for t, i in zip(temperatures, range(len(temperatures))):
            for m, j in zip(mass_NaCl_over_mass_water, range(len(mass_NaCl_over_mass_water))):
                x.append(t)
                y.append(m)
                d_table = densities.ravel()[count]
                d_func = a*t**2 + b*t + c*m**2 + d*m + e
                ztable.append(d_table)
                zfunc.append(d_func)
                if abs(d_table - d_func) > max_diff:
                    max_diff = abs(d_table - d_func)
                count += 1

        print "Density = %0.4g r^2 + %0.4g r + %0.4e t^2 + %0.4e t + %0.4f" % (c, d, a, b, e)
        print "Max diff = %0.4g" % max_diff
        print a, b, c, d, e
    elif choice.get() == "sucrose":
        molalities = numpy.array([0.015, 0.030, 0.060, 0.090, 0.122, 0.154, 0.186, 0.220, 0.254, 0.289, 0.325, 0.398, 0.476, 0.556, 0.641, 0.730, 0.824, 0.923, 1.026, 1.136, 1.252, 1.375, 1.505, 1.643, 1.791, 1.948, 2.116, 2.295, 2.489, 2.697, 2.921, 4.382, 6.817, 11.686])
        densities = numpy.array([1.0002, 1.0021, 1.0060, 1.0099, 1.0139, 1.0178, 1.0218, 1.0259, 1.0299, 1.0340, 1.0381, 1.0465, 1.0549, 1.0635, 1.0722, 1.0810, 1.0899, 1.0990, 1.1082, 1.1175, 1.1270, 1.1366, 1.1464, 1.1562, 1.1663, 1.1765, 1.1868, 1.1972, 1.2079, 1.2186, 1.2295, 1.2864, 1.3472, 1.4117])
        mass_sucrose_over_mass_water = molalities * 342.2965 / 1000.0

        def error_function((c, d, e)):
            sum = 0.0
            for m, j in zip(mass_sucrose_over_mass_water, range(len(mass_sucrose_over_mass_water))):
                sum += ((c*m**2 + d*m + e) - densities[j])**2
            return sum

        c, d, e = scipy.optimize.fmin(error_function, (3, 4, 1), full_output=1, xtol=1e-29, ftol=1e-29, maxiter=10000, maxfun=10000)[0]

        y = []
        ztable = []
        zfunc = []
        count = 0
        max_diff = 0.0

        for m, j in zip(mass_sucrose_over_mass_water, range(len(mass_sucrose_over_mass_water))):
            y.append(m)
            d_table = densities.ravel()[count]
            d_func = c*m**2 + d*m + e
            ztable.append(d_table)
            zfunc.append(d_func)
            if abs(d_table - d_func) > max_diff:
                max_diff = abs(d_table - d_func)
            count += 1

        print "Density = %0.4g r^2 + %0.4g r + %0.4f" % (c, d, e)
        print "Max diff = %0.4g" % max_diff
        return c, d, e


main_frame = Tkinter.Frame(root)
main_frame.pack()
w = Tkinter.Label(main_frame, text="Solute:")
w.pack(side=Tkinter.LEFT)
o = Tkinter.OptionMenu(main_frame, choice, "select...", "NaCl", "sucrose", command=o_changed)
o.pack(side=Tkinter.RIGHT)

solute_mass_frame = Tkinter.Frame(root)
solute_mass_frame.pack()

solute_mass_label = Tkinter.Label(solute_mass_frame, text="Solute mass:")
solute_mass_label.pack(side=Tkinter.LEFT)
solute_mass_entry = Tkinter.Entry(solute_mass_frame)
solute_mass_entry.pack(side=Tkinter.RIGHT)

solvent_mass_frame = Tkinter.Frame(root)
solvent_mass_frame.pack()

solvent_mass_label = Tkinter.Label(solvent_mass_frame, text="Solvent mass:")
solvent_mass_label.pack(side=Tkinter.LEFT)
solvent_mass_entry = Tkinter.Entry(solvent_mass_frame)
solvent_mass_entry.pack(side=Tkinter.RIGHT)

b = Tkinter.Button(root, text="pressme", command=callback)
b.pack()


root.mainloop()
