import sys
from tcr_complex import TCRComplex 

def main():
    
    tcrfile = ""
    pepmhcfile = ""
    mhc_type = 0
    xshift = 0.0
    crossshift = 0.0
    angleshift = 0.0
    shift_set = False

    
    argc = len(sys.argv)

    if argc < 2:
        print(
            "usage:\n"
            "tcr_docking_angle tcr_complex_file\n"
            "OR\n"
            "tcr_docking_angle tcr_file mhcpep_file\n"
            "OR\n"
            "tcr_docking_angle tcr_complex_file mhc_type\n"
            "OR\n"
            "tcr_docking_angle tcr_file mhcpep_file mhc_type\n"
            "MHC Type: 0 = Class I; 1 = Class II; 2 = CD1d, CD1b; "
            "3 = MR1; 4 = CD1d (4EI5); 5 = CD1c"
        )
        sys.exit(0)
    
    elif argc == 2:
        tcrfile = sys.argv[1]

    elif argc == 3:
        tcrfile = sys.argv[1]
        pepmhcfile = sys.argv[2]
        if len(pepmhcfile) == 1 and pepmhcfile.isdigit():
            mhc_type = int(pepmhcfile)
            pepmhcfile = ""

    elif argc == 4:
        tcrfile = sys.argv[1]
        pepmhcfile = sys.argv[2]
        mhc_type = int(sys.argv[3])

    elif argc == 6:
        tcrfile = sys.argv[1]
        mhc_type = int(sys.argv[2])
        xshift = float(sys.argv[3])
        crossshift = float(sys.argv[4])
        angleshift = float(sys.argv[5])
        shift_set = True

    elif argc == 7:
        tcrfile = sys.argv[1]
        pepmhcfile = sys.argv[2]
        mhc_type = int(sys.argv[3])
        xshift = float(sys.argv[4])
        crossshift = float(sys.argv[5])
        angleshift = float(sys.argv[6])
        shift_set = True

    # TCRComplex instance
    comp = TCRComplex(tcrfile, pepmhcfile, mhc_type)

    
    if shift_set:
        comp.set_x_shift(xshift)
        comp.set_cross_shift(crossshift)
        comp.set_tilt_shift(angleshift)

    
    comp.calc_docking_angle()

if __name__ == "__main__":
    main()
