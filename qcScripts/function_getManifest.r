
# function to download and load the manifest file

getManifest <- function(refDir = "."){
        if(!file.exists(file.path(refDir,"MouseMethylation-12v1-0_A2.csv"))){
            print("Downloading manifest...")
            manifest="https://emea.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/mouse-methylation/infinium-mouse-methylation-manifest-file-csv.zip"
            destfile <- "manifest.zip" 
            download.file(manifest, file.path(refDir,destfile))
            unzip(file.path(refDir,destfile), exdir = refDir)
            unlink(file.path(refDir,destfile))
            print("Manifest successfully downloaded")
        }

    print("Loading manifest...")
    return(data.table::fread(file.path(refDir,"MouseMethylation-12v1-0_A2.csv"), skip=7, fill=TRUE, data.table=F))

}
