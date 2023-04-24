from typing import List

from dnbc4tools.tools.utils import change_path, judgeFilexits
from dnbc4tools.__init__ import __root_dir__


def change_tag(inbam: str, outbam: str, tag1: str, tag2: str):
    """
    Changes the contents of two tags in a BAM file.

    Args:
        inbam: Input BAM file name/path.
        outbam: Output BAM file name/path.
        tag1: Name of first tag to be changed.
        tag2: Name of second tag to be changed.

    Returns:
        None
    """
    import pysam
    # Open the input BAM file and the output BAM file
    with pysam.AlignmentFile(inbam, "rb") as samfile, pysam.AlignmentFile(outbam, "wb", header=dict(samfile.header)) as outsam:
        
        # Loop over each read in the input BAM file
        for sam in samfile:
            try:
                
                # Check if the read has the specified tag
                if sam.has_tag('UB'):
                    Rtag1 = sam.get_tag(tag1)
                    Rtag2 = sam.get_tag(tag2)
                    
                    # Swap the contents of the two tags
                    sam.set_tag(tag1, Rtag2)
                    sam.set_tag(tag2, Rtag1)
                    
                    # Write the updated read to the output BAM file
                    outsam.write(sam)
            except KeyError:
                continue


class ChangeTag:
    def __init__(self, args):
        self.inbam: str = args.inbam
        self.outbam: str = args.outbam
        self.tag: str = args.tag

    def run(self):
        """
        Runs the Changetag command.

        Args:
            None

        Returns:
            None
        """
        # Check if input BAM file exists
        judgeFilexits(self.inbam)
        
        # Change the current working directory to the root directory of the package
        change_path()
        print("\033[0;32;40mStart Analysis\033[0m")
        
        # Split the specified tag into two parts
        tagsplit: List[str] = self.tag.split(',')
        print("Using \"%s\" for Analysis" % self.tag)
        if len(tagsplit) == 2:
            
            # Call the change_tag function to swap the contents of the two tags
            change_tag(self.inbam, self.outbam, tagsplit[0], tagsplit[1])
        else:
            raise Exception('The tag format or input is wrong.')
        print("\033[0;32;40mComplete\033[0m")


def changetag(args):
    """
    Runs the ChangeTag command.

    Args:
        args: Command line arguments for ChangeTag command.

    Returns:
        None
    """
    ChangeTag(args).run()


def helpInfo_changetag(parser):
    """
    Adds command line arguments to the parser for the ChangeTag command.

    Args:
        parser: The argparse.ArgumentParser object.

    Returns:
        None
    """
    parser.add_argument(
        '--inbam',
        metavar='FILE',
        help='BAM file that is used as input.',
        type=str,
        required=True
    )
    parser.add_argument(
        '--outbam',
        metavar='FILE',
        help='BAM file that is generated as output.',
        type=str,
        required=True
    )
    parser.add_argument(
        '--tag',
        metavar='TAGS',
        default='CB,DB',
        help='Tags that need to exchange content, [default: CB,DB]', 
        )
    return parser
