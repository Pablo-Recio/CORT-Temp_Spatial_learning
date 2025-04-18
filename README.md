# CORT-Temp_SpaL

#### Research paper

This repository contains the final code, data, and figures used in the following manuscript:

Recio et al. 2024. Prenatal conditions do not affect brain physiology and learning in a lizard

Supplementary Materials: The supplementary materials associated with this paper is integrated within the ms.docx or ms.qmd files. 

#### How to use this repository?

Users can download a zip file of the entire repository by clicking on the green code tab at the top of the page and then clicking Download ZIP. Users who already have a GitHub account can fork the repository.

The key file in this repository is the 📄 ms.qmd. This file can be rendered in R with Quarto to reproduce the entire paper. Code chunks within the file provide the code used to reproduce figures and analyses along with supporting statements within the text. Note that inline code chunks use specific objects which are then rendered.

The 📄 ms.qmd file makes use of files within a number of folders that are identified in the code chunks. There are a number of important folders in the repository.

  📂 data folder contains all the raw data used in files. Note that there are different files, but the file metadata.xlsx contains all the information about each of the files. In output/database_clean there are the main files employed in the analyses and figures; names of the columns are defined in the metadata.xls file on data folder.   
  📂 output/figs/ Folder contains all the figures for the paper that are read and included in the paper.   
  📂 R folder contains three files, two of them used to clean and process data to prepare it for use in the ms.qmd file. Note that readers do not need to open and run these files, but they are simply here to document the workflow and code used to clean up data to be used. These include:  
        📄 1_data_process.R, which is used to produce the files in output/database_clean from the raw data in data folder.  
        📄 func.R, which contains all the functions called later in the ms.qmd file.   
  📂 bib The bib folder contains:  
        📄 refs.bib the references;  
        📄 animal-cognition.csl the journal formatting style file;  
        📄 template.docx a template docx file to format the resulting rendered files.  
  📂 ANIMCOG This folder contains all the material submitted to the journal, plus all the material from the reviews.  

#### Reporting Issues or Asking Questions
If anything is unclear or you require further detail please do not hesitate to lodge an issue.
