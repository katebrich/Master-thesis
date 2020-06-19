import sys
import getopt

try:
    opts, args = getopt.getopt(sys.argv[1:], 'd:f:')
except getopt.GetoptError as err:
    # print help information and exit:
    print(err) # will print something like "option -a not recognized"
    #todo print help
    sys.exit(2)
for opt, arg in opts:
    if opt in ("-h", "--help"):
        # todo print help
        sys.exit()
    elif opt in ("-d", "--directory"):
        dir = arg
    elif opt in ("-f", "--file"):
        file = arg
    #todo else unknown option

