#
# Simple cross-section for O16 using mostly defaults and 
# FIXME: wrong title and missing files in repo for custom wt 
#
CWD=$PWD

Isotopes=(Al27 Ar36)
#(Al27, Ar36, Ar38, Ar40, C12, C13, Cd106, Cd108, Cd110, Cd111, Cd112, Cd114, Cd116, Cr52, Cu63, Cu65, Fe56, H1, He3, Mg24, Mg25, Mn55, N14, N15, O16, O17, Si28, Si29, Si30, Ti48, Zn64, Zn66, Zn67, Zn68, ZN70)
            
for i in Fe56 Al27 C12 Cd106 Cd108 Cd110 Cd111 Cd112 Cd113 Cd114 Cd116 Cr52 H1 He3 N14 Ni58 O16 Si28 Ar36 Ar38 Ar40 Cu63 Cu65 Mg26 Mn55 Ti48 Zn64 Zn66 Zn67 Zn68 Zn70
do
   echo "adding permission for $i.sh"
   chmod +x $i.sh 
done

cd "$CWD" || exit



