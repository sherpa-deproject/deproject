#
# This assumes fit_m87.py has been run and that no significant
# changes have been made.
#
# It is run as
#    %run -i fit_m87_tied.py
#

print(dep.get_shells())

dep.tie_par('xsmekal.kt', 17, 18)
dep.tie_par('xsmekal.abundanc', 17, 18)

print(dep.get_shells()[0:5])

dep.fit()

tied = dep.get_fit_results()
print(len(tied))

print(tied['xsmekal.kT'][16:])
print(tied['xsmekal.Abundanc'][16:])

dep.fit_plot('xsmekal.kt', units='pc')
rmid = (onion['rlo_phys'] + onion['rhi_phys']) / 2
plt.scatter(rmid.to(u.pc), onion['xsmekal.kT'], c='k')
plt.savefig('m87_temperature_tied_comparison.png')

dep.untie_par('xsmekal.kt', 18)
dep.untie_par('xsmekal.abundanc', 18)
print(dep.get_shells()[0:4])
