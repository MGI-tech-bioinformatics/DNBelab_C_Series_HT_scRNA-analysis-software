<div align="right">

[🏠 Home](../README.md)

</div>

<br>

# DNBelab C Series Demo Datasets

<br>

Demo datasets are hosted on CNGB (China National GeneBank). Only mouse samples data is provided.

<br>

---

<br>

<table style="width:100%; border-collapse: collapse; margin: 1.5em 0;">
  <thead style="background-color: #f2f2f2; border-bottom: 2px solid #ddd;">
    <tr>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Analysis Type</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Sample Information</th>
      <th style="padding: 12px 15px; border: 1px solid #ddd; text-align: left;">Project Link</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>scRNA-seq v3</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Mouse sample with separate cDNA and Oligo libraries.</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><a href="https://db.cngb.org/data_resources/project/CNP0005575/" target="_blank">CNP0005575</a></td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>scVDJ-seq</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Mouse spleen tissue, including 5' RNA, TCR, and BCR data.</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><a href="https://db.cngb.org/data_resources/project/CNP0006116/" target="_blank">CNP0006116</a></td>
    </tr>
    <tr>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><b>scATAC-seq</b></td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;">Mouse brain tissue.</td>
      <td style="padding: 12px 15px; border: 1px solid #ddd;"><a href="https://db.cngb.org/data_resources/project/CNP0004369" target="_blank">CNP0004369</a></td>
    </tr>
  </tbody>
</table>

---

## Notes on Data Usage

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
💡 <strong>scVDJ-seq Tip</strong>: If you don't need to analyze the 5' RNA data, you can directly download the `singlecell.csv` file from the 5' RNA data directory on the FTP server to use as input for the VDJ pipeline.
</div>

<div style="background-color: #e7f3fe; border-left: 6px solid #2196F3; padding: 15px; margin: 1.5em 0; border-radius: 4px;">
💡 <strong>General Tip</strong>: For all datasets, you can click the project link, navigate to the "Sample" or "Experiment" tabs, and find the FTP link to download the raw data for analysis.
</div>

<div align="center">
  <img src="./images/scrna_v3_datasets.jpg" alt="scRNA-seq v3 Sample Details" width="800">
  <p><i>Example: Finding the FTP download link in the CNGB project page.</i></p>
</div>