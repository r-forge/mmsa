#######################################################################
# BiostringsTools - Interfaces to several sequence alignment and 
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

BiostringsTools_Software_Wizard <- function(rdp=FALSE, clustal=FALSE, kalign=FALSE, 
  mafft=FALSE, muscle = FALSE, blast=FALSE, blast16S=FALSE, 
  boxshade=FALSE) {
  
  os <- toupper(Sys.info()["sysname"])
  ### note: DARWIN is OSX
  
  cat("BiostringsTools Software Installation Wizard for", os, "\n\n")
  
  ### default install dir
  dir <- path.expand("~/BiostringsTools/")
  if(!file.exists(dir)) dir.create(dir)
  on.exit({ setwd(getwd()) })
  setwd(dir)
  
  ### do all (not blast16S)?
  if(!(any(c(rdp, clustal, kalign, mafft, muscle,
    blast, blast16S, boxshade)))) {
    rdp <- TRUE; clustal<- TRUE; kalign <- TRUE; mafft <- TRUE; muscle <- TRUE;
    blast <- TRUE; blast16S <- FALSE; boxshade <- TRUE
  }
  
  if(rdp) {
    if(!file.exists("rdp_classifier_2.5")) {
      cat("Installing RDP\n")
      ### any OS
      download.file("http://downloads.sourceforge.net/project/rdp-classifier/rdp-classifier/rdp_classifier_2.5.zip?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Frdp-classifier%2Ffiles%2Frdp-classifier%2F&use_mirror=iweb",
        "rdp_classifier_2.5.zip", mode='wb')
      unzip("rdp_classifier_2.5.zip")
      Sys.setenv(RDP_JAR_PATH=file.path(dir,"rdp_classifier_2.5","rdp_classifier-2.5.jar"))
    }else{
      cat("RDP ... installed.\n")
    }
  }
  
  if(clustal) {
    if(any(Sys.which(c("clustalw", "clustalw2"))!="")) cat("clustalw ... installed.\n")
    else if(os=="WINDOWS") {
      cat("Installing clustal\n")
      download.file("http://www.clustal.org/download/current/clustalw-2.1-win.msi",
        "clustalw-2.1.msi", mode='wb')
      system("msiexec /i clustalw-2.1.msi")
    } else if(os=="DARWIN") {
      download.file("http://www.clustal.org/download/current/clustalw-2.1-macosx.dmg",file.path(dir,"clustalw-2.1-macosx.dmg"),mode='wb')
      system(paste("hdiutil mount",file.path(dir,"clustalw-2.1-macosx.dmg")))
      system("cp -R /Volumes/clustalw-2.1-macosx /Applications")
      system("hdiutil unmount /Volumes/clustalw-2.1-macosx")
    } else cat("Please install package 'clustalw' manually (in package manager)!\n")
  }
  
  if(kalign) { 
    if(any(Sys.which("kalign")!="")) cat("kalign ... installed.\n")
    else if(os=="Windows") {
      cat("Installing kalign\n")
      download.file("http://msa.sbc.su.se/downloads/kalign/current.tar.gz",
        "kalign.tar.gz", mode='wb')
      untar("kalign.tar.gz")
      cat("Please follow instructions in the file ",file.path(dir,"kalign","README"))
    } else if(os=="DARWIN") {
      cat("Installing kalign\n")
      download.file("http://msa.sbc.su.se/downloads/kalign/current.tar.gz",file.path(dir,"kalign.tar.gz"),mode='wb')
      untar(file.path(dir,"kalign.tar.gz"),exdir=file.path(dir,"kalign"))
      setwd("kalign")
      system("./configure")
      system("make")
      system("osascript -e 'do shell script \"sudo make install\" with administrator privileges'")
      #system(paste("sudo make install"))
      setwd("..")	
    } else
      cat("Please install 'kalign' manually (in package manager)!\n")
  }
  
  if(mafft) {
    if(any(Sys.which("mafft")!="")) cat("mafft ... installed.\n")
    else if(os=="Windows") {	
      cat("Please install package 'mafft' manually!\n")
    } else if(os=="DARWIN") {
      download.file("http://mafft.cbrc.jp/alignment/software/mafft-7.050-signed.pkg",file.path(dir,"mafft-7.050-signed.pkg"),mode='wb')
      cmd <- paste("installer -pkg",file.path(dir,"mafft-7.050-signed.pkg"),"-target /")
      #system("osascript -e 'do shell script \"installer -pkg\" & 
      system(cmd)
    }else
      cat("Please install 'mafft' manually (in package manager)!\n")
  }
  
  if(blast) {
    if(any(Sys.which("blastn")!="")) cat("BLAST ... installed.\n")
    else if(os=="WINDOWS") {
      cat("Installing BLAST\n")
      download.file("ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.28+-win64.exe",
        "ncbi-blast-2.28+-win64.exe", mode='wb')
      system("ncbi-blast-2.28+-win64.exe")
    } else if(os=="DARWIN") {
      download.file("ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.2.28+.dmg",file.path(dir,"ncbi-blast-2.2.28+.dmg"),mode='wb')
      system("sudo hdiutil mount ncbi-blast-2.2.28+.dmg")
      system("sudo cp -R /Volumes/ncbi-blast-2.2.28+/ /Applications")
      system("sudo hdiutil unmount ncbi-blast-2.2.28+")    
    } else cat("Please install package 'ncbi-blast+' manually!\n")
  }
  
  if(blast16S) {
    if(file.exists("16SMicrobialDB")) cat("16SMicrobialDB ... installed.\n")
      else {
        cat("Installing the 16S rRNA database for BLAST\n")
        download.file("ftp://ftp.ncbi.nlm.nih.gov/blast/db/16SMicrobial.tar.gz",
          "16SMicrobial.tar.gz", mode='wb')
        untar("16SMicrobial.tar.gz", exdir="16SMicrobialDB")
      }
  }
  
  if(boxshade) {
    if(any(Sys.which("boxshade")!="")) cat("boxshade ... installed.\n")
    else cat("Please install package 'boxshade' manually!\n")
  }
  
  if(muscle) {
    if(any(Sys.which("muscle")!="")) cat("MUSCLE ... installed.\n")
    else cat("Please install package 'muscle' manually!\n")
  }
}
