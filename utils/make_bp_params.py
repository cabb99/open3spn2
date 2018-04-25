#!/usr/bin/env python

""" 
Scripts for generating base step parameter files for reading by 
X3DNA.  Reads a sequence file formatted as follows:

<Number of Base Pairs>
<Sequence>

e.g. 
20
TTGCCAACGTCCGACTGAAA
"""

import sys
from BaseSteps import *

# Defining a few constant values
baseDict = {"A":0,"T":1,"C":2,"G":3}
FOUR = 4

# Read the sequence file

def get_command_line_args(args):
    if len(args) != 2:
        print "Usage: %s <sequence file>" % args[0]
        sys.exit(1)
    return args[1]

def read_sequence_file(fileName):
    try:
        file = open(fileName,"r")
        try:
            N = int(file.readline())
            sequence = file.readline().strip()
            if (N != len(sequence)):
                print "Inconsistent number of bases in the sequence file"
                sys.exit(2)
            else:
                # Check sequence file for spurious values
                for letter in sequence:
                    if letter not in ["A","T","C","G"]:
                        print "Unrecognized letter in sequence file (%s)" % letter
                        sys.exit(2)
                return sequence

        except ValueError:
            print "Number of bases is not a number"
            sys.exit(2)

    except IOError:
        print "Count not open %s" % filename
        sys.exit(2)


def write_bp_parameters(sequence):
    outFile = "bp_step.par"
    counter = 0
    last = 100
    avg_twist = 0
    N = len(sequence)

    # Write the header information
    file = open(outFile,"w")
    file.write("%4d # base-pairs\n" % N)
    file.write("   0 # ***local base-pair & step parameters***\n")
    file.write("#        Shear    Stretch   Stagger   Buckle   Prop-Tw   Opening     Shift     Slide     Rise      Tilt      Roll      Twist\n")
    
    # Write the base step parameters
    for i in range(N):
        bpID = baseDict[sequence[i]]
        mybp = bp[bpID]
        file.write("%s-%s%11.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf" \
            % (mybp.stea,mybp.steb,mybp.shear,mybp.stretch,mybp.stagger,mybp.buckle,mybp.propeller,mybp.opening))
        if i == 0:
            file.write("%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf\n" \
                % (0.0,0.0,0.0,0.0,0.0,0.0))
        else:
            bsID = baseDict[sequence[i]]
            mybs = step[last*FOUR+bsID]
            file.write("%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf%10.3lf\n" \
                % (mybs.shift,mybs.slide,mybs.rise,mybs.tilt, mybs.roll, mybs.twist))
            counter += 1
            avg_twist += mybs.twist
        last = bpID

    avg_twist /= float(counter)
    print "The average twist in this sequence is %lf\n" % avg_twist
    file.close()

def main():
    sequenceFile = get_command_line_args(sys.argv)
    sequence = read_sequence_file(sequenceFile)
    write_bp_parameters(sequence)

if __name__ == "__main__":
    main()






