import os
import shutil
from typing import Any
from dnbc4tools.tools.utils import str_mkdir, judgeFilexits, change_path, read_json
from dnbc4tools.__init__ import __root_dir__

class Report:
    def __init__(self, args: Any) -> None:
        self.name: str = args.name
        self.outdir: str = os.path.join(args.outdir, self.name)
        self.genomeDir: str = args.genomeDir

    def run(self) -> None:
        judgeFilexits(self.genomeDir)
        str_mkdir('%s/04.report/div'%self.outdir)
        str_mkdir('%s/04.report/base64'%self.outdir)
        # str_mkdir('%s/04.report/table'%self.outdir)
        str_mkdir('%s/log'%self.outdir)
        str_mkdir('%s/output'%self.outdir)
        change_path()

        genomeDir: str = os.path.abspath(self.genomeDir)
        indexConfig: Any = read_json('%s/ref.json'%genomeDir)
        species: str = indexConfig['species']

        print('\nSummarizing results')
        from dnbc4tools.tools._plotly import plot_barcode_atac_frag,plot_jaccard_atac_frag,atac_saturation,plot_cluster,plot_cluster_uniqueFrags
        judgeFilexits(
            '%s/02.decon/%s.barcodeCount.tsv'%(self.outdir,self.name),
            '%s/02.decon/%s.CorrelationBarcodes.tsv.gz'%(self.outdir,self.name),
            '%s/02.decon/%s.d2cCutoff.tsv'%(self.outdir,self.name),
            '%s/02.decon/%s.sequenceSaturation.tsv'%(self.outdir,self.name),
            '%s/03.analysis/cluster_cell.stat'%self.outdir
        )
        plot_barcode_atac_frag(
            '%s/02.decon/%s.barcodeCount.tsv'%(self.outdir,self.name),
            '%s/02.decon/%s.d2cCutoff.tsv'%(self.outdir,self.name),
            '%s/04.report/div'%self.outdir
            )
        plot_jaccard_atac_frag(
            '%s/02.decon/%s.CorrelationBarcodes.tsv.gz'%(self.outdir,self.name),
            '%s/02.decon/%s.d2cCutoff.tsv'%(self.outdir,self.name),
            '%s/04.report/div'%self.outdir
        )
        atac_saturation(
            '%s/02.decon/%s.sequenceSaturation.tsv'%(self.outdir,self.name),
            '%s/04.report/div'%self.outdir
        )
        plot_cluster(
            '%s/03.analysis/cluster_cell.stat'%self.outdir,
            '%s/04.report/div'%self.outdir
        )
        plot_cluster_uniqueFrags(
            '%s/03.analysis/cluster_cell.stat'%self.outdir,
            '%s/04.report/div'%self.outdir
        )
        from dnbc4tools.tools.utils import png_to_base64
        judgeFilexits(
            '%s/03.analysis/images/DropBeadsnum.png'%self.outdir,
            '%s/03.analysis/images/QC.png'%self.outdir,
            '%s/03.analysis/images/TSS.png'%self.outdir,
            '%s/03.analysis/images/InterSize.png'%self.outdir,
        )
        pictures = {
            '%s/03.analysis/images/DropBeadsnum.png'%self.outdir:'plot3_DropBeadsnum',
            '%s/03.analysis/images/QC.png'%self.outdir:'plot4_QC',
            '%s/03.analysis/images/InterSize.png'%self.outdir:'plot5_InterSize',
            '%s/03.analysis/images/TSS.png'%self.outdir:'plot6_TSS'
            }
        for infile,outfile in pictures.items():
            png_to_base64(infile,outfile,'%s/04.report/base64'%self.outdir)

        from dnbc4tools.atac.src.generate_report import write_param_to_template
        report_file, metrics_summary_df = write_param_to_template(
            '%s/template/template_scATAC.html'%__root_dir__,
            self.name,
            self.outdir,
            self.genomeDir
            )
        htmlFile = open('%s/04.report/%s_scATAC_report.html'%(self.outdir,self.name),'w')
        htmlFile.write(report_file)
        htmlFile.close()
        metrics_summary_df.to_csv('%s/04.report/metrics_summary.xls'%self.outdir,sep='\t',index=None)

        if os.path.exists('%s/04.report/metrics_summary.xls'%self.outdir):
            shutil.copy("%s/04.report/metrics_summary.xls"%self.outdir,'%s/output'%self.outdir)
        if os.path.exists('%s/04.report/%s_scATAC_report.html'%(self.outdir,self.name)):
            shutil.copy("%s/04.report/%s_scATAC_report.html"%(self.outdir,self.name),'%s/output'%self.outdir)
        if os.path.exists('%s/02.decon/%s.fragments.tsv.gz'%(self.outdir,self.name)):
            shutil.copy('%s/02.decon/%s.fragments.tsv.gz'%(self.outdir,self.name),'%s/output'%self.outdir)
        if os.path.exists('%s/02.decon/%s.fragments.tsv.gz.tbi'%(self.outdir,self.name)):
            shutil.copy('%s/02.decon/%s.fragments.tsv.gz.tbi'%(self.outdir,self.name),'%s/output'%self.outdir)
        if os.path.exists('%s/02.decon/%s.Metadata.tsv'%(self.outdir,self.name)):
            shutil.copy('%s/02.decon/%s.Metadata.tsv'%(self.outdir,self.name),'%s/output'%self.outdir)
        if os.path.exists("%s/03.analysis/peak"%self.outdir):
            if os.path.exists('%s/output/Peak_matrix'%self.outdir):
                shutil.rmtree('%s/output/Peak_matrix'%self.outdir)
            shutil.copytree("%s/03.analysis/peak"%self.outdir,'%s/output/Peak_matrix'%self.outdir,dirs_exist_ok=True)

def report(args):
    Report(args).run()


def helpInfo_report(parser):
    parser.add_argument(
        '--name', 
        metavar='NAME',
        help='Sample name.', 
        type=str,
        required=True
        )
    parser.add_argument(
        '--outdir', 
        metavar='PATH',
        help='Output diretory, [default: current directory].', 
        default=os.getcwd()
        )
    parser.add_argument(
        '--genomeDir',
        type=str, 
        metavar='PATH',
        help='Path of folder containing reference database.',
        required=True
        )
    return parser