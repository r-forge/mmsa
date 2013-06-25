#######################################################################
# BioTools - Interfaces to several sequence alignment and 
# classification tools
# Copyright (C) 2012 Michael Hahsler and Anurag Nagar
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

BioTools_Software_Wizard <- function(RDP=FALSE, clustal=FALSE, kalign=FALSE, 
                                     MAFFT=FALSE, BLAST=FALSE, BLAST16S=FALSE, 
                                     boxshade=FALSE) {

  os <- toupper(Sys.info()["sysname"])
  ### note: DARWIN is OSX
  
  cat("BioTools Software Installation Wizard for", os, "\n")
  
  ### default install dir
  dir <- path.expand("~/BioTools/")
  if(!file.exists(dir)) dir.create(dir)
  on.exit({ setwd(getwd()) })
  setwd(dir)
 
  ### do all (not BLAST16S)?
  if(!(any(c(RDP, clustal, kalign, MAFFT, 
           BLAST, BLAST16S, boxshade)))) {
    RDP <- TRUE; clustal<- TRUE; kalign <- TRUE; MAFFT <- TRUE; 
            BLAST <- TRUE; BLAST16S <- FALSE; boxshade <- TRUE
  }
  
  if(RDP) {
    if(!file.exists("rdp_classifier_2.5")) {
      cat("Installing RDP\n")
      ### any OS
      download.file("http://downloads.sourceforge.net/project/rdp-classifier/rdp-classifier/rdp_classifier_2.5.zip?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Frdp-classifier%2Ffiles%2Frdp-classifier%2F&use_mirror=iweb",
                    "rdp_classifier_2.5.zip", mode='wb')
      unzip("rdp_classifier_2.5.zip")
    }
  }
  
  if(clustal) {
    ### FIXME: check if it is already installed
    if(os=="WINDOWS") {
      cat("Installing clustal\n")
      download.file("http://www.clustal.org/download/current/clustalw-2.1-win.msi",
                    "clustalw-2.1.msi", mode='wb')
      system("msiexec /i clustalw-2.1.msi")
    } else if(os=="DARWIN") {
      cat("Please install clustal manually from http://www.clustal.org/download/current/\n")
      #download.file("http://www.clustal.org/download/current/clustalw-2.1-macosx.dmg",
      #              "clustalw-2.1-macosx.dmg",mode='wb')
      #system("hdiutil mount clustalw-2.1-macosx.dmg")
      #system("sudo cp -R /Volumes/clustalw-2.1-macosx /Applications")
      #system("hdiutil unmount clustalw-2.1-macosx")
    } else cat("Please install package 'clustal' manually (in package manager)!\n")
  }
  
  if(kalign) { 
    ### FIXME: check if it is already installed
    if(os=="Windows") {
      cat("Please install 'kalign' manually!\n")
    } else if(os=="DARWIN") {
      cat("Installing kalign\n")
      download.file("http://msa.sbc.su.se/downloads/kalign/current.tar.gz",
                    "kalign.tar.gz", mode='wb')
      untar("kalign.tar.gz", exdir="kalign")
      system("kalign/configure")
      system("kalign/make")
      system(paste("sudo make","install"))
    }
    cat("Please install 'kalign' manually (in package manager)!\n")
  }
  
  if(MAFFT) {
    cat("Please install package 'mafft' manually!\n")
  }
  
  if(BLAST) {
    ### FIXME: check if it is already installed
    if(os=="WINDOWS") {
      cat("Installing BLAST\n")
      download.file("ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.28+-win64.exe",
                    "ncbi-blast-2.28+-win64.exe", mode='wb')
      System("ncbi-blast-2.28+-win64.exe")
    }else if(os=="DARWIN") {
      cat("Install BLAST manually from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/\n")
      #download.file("ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.28+.dmg",file.path(downloadDir,"ncbi-blast-2.2.28+.dmg"),mode='wb')
      #system("hdiutil mount ncbi-blast-2.2.28+.dmg")
      #system("sudo cp -R /Volumes/ncbi-blast-2.2.28+/ /Applications")
      #system("hdiutil unmount ncbi-blast-2.2.28+")
    } else cat("Please install package 'blast+' manually!\n")
  }
  
  if(BLAST16S) {
    if(file.exists("16SMicrobialDB"))
    cat("Installing the 16S rRNA database for BLAST\n")
    download.file("ftp://ftp.ncbi.nlm.nih.gov/blast/db/16SMicrobial.tar.gz",
                  "16SMicrobial.tar.gz", mode='wb')
    untar("16SMicrobial.tar.gz", exdir="16SMicrobialDB")
  }
  
  if(boxshade) {
    cat("Please install package 'boxshade' manually!\n")
  }
}
