import argparse
from scripts.StrainIQ_builder import *
from scripts.StrainIQ_identifier import *
from scripts.StrainIQ_quantifier import *

def main():
    my_parser = argparse.ArgumentParser(description='List the content of a folder')
    my_parser.add_argument('-p',
                        metavar='Program',
                        type=str,
                        choices=['builder', 'identifier', 'quantifier'],
                        help='Program to run. Example: builder, identifier, quantifier',
                        required=True)
    my_parser.add_argument('-n',
                       metavar='n-size',
                       type=int,
                       help='Size of the n-gram',
                       default=21)
    my_parser.add_argument('-glist',
                       metavar='file with a list of reference genome',
                       type=str,
                       help='Tab delimited file with list of genomes to be included in the DSEM. Format: gid    filelocation')
    my_parser.add_argument('-ng',
                        metavar='total number of genomes in the dasem',
                        type=int,
                        help='The number of genomes in the DSEM. This parameteris needed if the DSEM is rebuilt using different sets for references.',
                        default=471)
    my_parser.add_argument('-dsem',
                        metavar='DSEM',
                        type=str,
                        help='The DSEM for this bodysite.')
    my_parser.add_argument('-sample',
                        metavar='sample',
                        type=str,
                        help='Input metagenomic sample. Convert the reverse reads in case of paired end sequencing before bombining the reads together')
    my_parser.add_argument('-prediction',
                        metavar='-prediction',
                        type=str,
                        help='The prediction file produced by the identifier')

    args = my_parser.parse_args()
    if args.p == "builder":
        if args.glist is None:
            my_parser.error("builder requires -glist")
            sys.exit()
        if not os.path.isfile(args.glist):
            print('The path specified does not exist')
            sys.exit()
        model(args.n, args.glist)
        print("builder")
    elif args.p == "identifier":
        #print(args.p)
        #print(args.n)
        #print(args.ng)
        if args.dsem is None:
             my_parser.error("identifier require -dsem")
        if args.sample is None:
            my_parser.error("identifier require -sample")
        #print(args.dsem)
        #print(args.sample)
        identifier(args.dsem, args.sample)
    elif args.p == "quantifier":
        if args.dsem is None:
            my_parser.error("identifier require -dsem")
        if args.sample is None:
            my_parser.error("identifier require -sample")
        if args.prediction is None:
            my_parser.error("identifier require -prediction")
        quantifier(args.dsem, args.sample,args.prediction)
        #parse_score_file(model, prediction) # generate identified genomes based on the cutoff

        print(args.p)



if __name__ == "__main__":
    main()

