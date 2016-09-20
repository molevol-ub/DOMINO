#!/bin/sh
## Script for unzipping and installing several binaries in DOMINO package
(
installer=$0
current_dir=$(dirname $installer)
install_dir=$(dirname $current_dir)

if [ $current_dir = "." ] 
then
	current_dir=`pwd`
	install_dir=$(dirname $current_dir)
fi

echo ""
echo "##############################################"
echo "	DOMINO installer"
echo "##############################################"
echo ""
echo "+ Working dir: "$current_dir
scripts_folder=$install_dir"/"scripts
mkdir $scripts_folder
## generate a bin folder
mkdir $install_dir"/"bin
mkdir $install_dir"/src"
qt_code=$install_dir"/src/Qt"
mkdir $qt_code
perl_code=$install_dir"/src/perl"
mkdir $perl_code
echo "+ Generate a folder for DOMINO scripts: "$scripts_folder
tmp_dir_name="tmp_dir"
tmp_dir_path=$install_dir"/"$tmp_dir_name
echo "+ Generate an intermediate folder: "$tmp_dir_path
mkdir $tmp_dir_path
echo "+ Copy, extract and compile necessary files..."
binaries_tar=$current_dir"/"third_party_binaries-src.tar-Mac.gz
cp $binaries_tar $tmp_dir_path
cd $tmp_dir_path
tar -zxvf "./"third_party_binaries-src.tar-Mac.gz
cd binaries
ls | while read folder_tar_gz
do
	echo "	- Extracting and compiling: "$folder_tar_gz
	cd $folder_tar_gz
	ls | while read tar_gz
	do
		tar -zxvf $tar_gz
		rm $tar_gz

		if [ "$folder_tar_gz" = "CAP3" ];
		then
			mkdir -p $scripts_folder"/"cap3"/"bin
			mv cap3.macosx.intel64"/"cap3 $scripts_folder"/"cap3"/"bin
			mv cap3.macosx.intel64"/"README $scripts_folder"/"cap3
			mv cap3.macosx.intel64"/"doc $scripts_folder"/"cap3
		fi
		if [ "$folder_tar_gz" = "bowtie" ];
		then
			mkdir -p $scripts_folder"/"bowtie2-2.2.9
			cd bowtie2-2.2.9
			make > DOMINO_Error.log			
			mv bowtie2* $scripts_folder"/"bowtie2-2.2.9
			mv AUTHORS $scripts_folder"/"bowtie2-2.2.9
			mv LICENSE $scripts_folder"/"bowtie2-2.2.9
			mv VERSION $scripts_folder"/"bowtie2-2.2.9
		fi

		if [ "$folder_tar_gz" = "mothur" ];
		then 
			mkdir -p $scripts_folder"/"MOTHUR_v1.32.0
			mv mothur"/"mothur $scripts_folder"/"MOTHUR_v1.32.0
			mv mothur"/"LICENSE $scripts_folder"/"MOTHUR_v1.32.0
		fi
		if [ "$folder_tar_gz" = "prinseq" ];
		then
			mkdir -p $scripts_folder"/"PRINSEQ-lite_0.20.4
			mv prinseq-lite-0.20.4"/"prinseq-lite.pl $scripts_folder"/"PRINSEQ-lite_0.20.4			
			mv prinseq-lite-0.20.4"/"README $scripts_folder"/"PRINSEQ-lite_0.20.4			
			mv prinseq-lite-0.20.4"/"COPYING $scripts_folder"/"PRINSEQ-lite_0.20.4
		fi

		if [ "$folder_tar_gz" = "MIRA" ]
		then
			mkdir -p $scripts_folder"/"mira_v4.0"/"bin
			mira_folder=mira_4.0.2_darwin13.1.0_x86_64_static
			mv $mira_folder"/"bin"/"mira $scripts_folder"/"mira_v4.0"/"bin
			mv $mira_folder"/"LICENCE $scripts_folder"/"mira_v4.0
			mv $mira_folder"/"README $scripts_folder"/"mira_v4.0
			mv $mira_folder"/"lib $scripts_folder"/"mira_v4.0
		fi			

		if [ "$folder_tar_gz" = "NGSQC" ]
		then
			mkdir -p $scripts_folder"/"NGSQCToolkit_v2.3.1
			mv NGSQCToolkit_v2.3.3"/"* $scripts_folder"/"NGSQCToolkit_v2.3.1
		fi

		if [ "$folder_tar_gz" = "NCBI" ]
		then
			mkdir -p $scripts_folder"/NCBI_BLAST"
			cd blast
			ls | while read files
			do
				cp $files $scripts_folder"/NCBI_BLAST"
			done
		fi
		
		if [ "$folder_tar_gz" = "samtools" ]
		then
			mkdir -p $scripts_folder"/"samtools-1.3.1
			cd samtools-1.3.1
 			make > DOMINO_Error.log
			mv COPYING $scripts_folder"/"samtools-1.3.1
			mv AUTHORS $scripts_folder"/"samtools-1.3.1
			mv samtools $scripts_folder"/"samtools-1.3.1
		fi			
	done
	cd $tmp_dir_path"/"binaries
done
rm -rf $tmp_dir_path
cd $current_dir
echo ""
echo "+ All necessary binaries have been extracted and compiled..."
echo ""
echo "##############################################"
echo "	Perl modules"
echo "##############################################"
echo ""
echo "+ Installing some perl modules:"
tar -zxvf tarGZ_modules_Perl_NGSQCToolkit.tar.gz
perl $current_dir"/"unzip_install_NGSQCtoolkit_modules.pl $current_dir"/"tarGZ_modules_Perl_NGSQCToolkit $install_dir"/"scripts"/"NGSQCToolkit_v2.3.1
tar -zxvf tarGZ_modules_Perl.tar.gz
perl $current_dir"/"unzip_install_perl_modules.pl $current_dir"/"tarGZ_modules_Perl $install_dir"/"scripts
cp -r $install_dir"/"scripts"/"NGSQCToolkit_v2.3.1"/"lib"/"Parallel $install_dir"/"scripts"/"lib
cp -r $scripts_folder"/"NGSQCToolkit_v2.3.1"/"lib"/"Parallel $scripts_folder"/"lib

## Deleting temporary files
cd $current_dir
rm -rf tarGZ_modules_Perl.tar.gz
rm -rf tarGZ_modules_Perl_NGSQCToolkit.tar.gz
rm -rf tarGZ_modules_Perl
rm -rf tarGZ_modules_Perl_NGSQCToolkit
rm unzip_install_perl_modules.pl
rm unzip_install_NGSQCtoolkit_modules.pl
rm -rf binaries_tar

echo ""
echo "+ Perl modules have been extracted into different folders..."
echo ""
echo ""
echo "##############################################"
echo "	Finishing installation"
echo "##############################################"
echo ""
echo "+ Copying DOMINO scripts"
cp DM_*pl $scripts_folder
cp DOMINO.pm $scripts_folder"/"lib
cp LICENSE.txt $scripts_folder
cp DM_*pl $perl_code
mv db_default.tar.gz $scripts_folder
cd $scripts_folder
tar -zxvf db_default.tar.gz
rm -r db_default.tar.gz
cd $current_dir

#copy necesary libs and fonts...
mv $current_dir"/"fonts.tar.gz $install_dir"/"bin
cd $install_dir"/"bin
tar -zxvf fonts.tar.gz
rm fonts.tar.gz
cd $current_dir
mv $current_dir"/Qt-c++/"DOMINO.app $install_dir"/"bin
mv $current_dir"/Qt-c++/"* $qt_code
mv $current_dir"/"starter_Mac.command $install_dir"/"bin
mv $install_dir"/"scripts $install_dir"/"bin
ln -s $install_dir"/"bin"/"DOMINO.app ~/Desktop/
ln -s $install_dir"/"bin"/"DOMINO.app $install_dir
mv uninstaller.command $install_dir"/"bin

echo "+ Copying necessary files..."

mkdir $install_dir"/"logo
mv $current_dir"/"*png $install_dir"/"logo
mv $current_dir"/"*ico $install_dir"/"logo
mv $current_dir"/"README.md $install_dir"/"
mv $current_dir"/"NEWS $install_dir"/"
mv $current_dir"/"TUTORIAL $install_dir"/"
mv $current_dir"/"VERSION $install_dir"/"
mv $current_dir"/"Change.log $install_dir"/"
mv $current_dir"/"LICENSE.txt $install_dir"/"
mv $current_dir"/"example $install_dir"/"
cd $install_dir

echo "+ Removing temporary files..."
cd $install_dir
rm -r $current_dir

#import variables
export LD_LIBRARY_PATH=$current_dir"/"lib
export QT_QPA_FONTDIR=$current_dir"/"fonts
chmod -R u+rwx $install_dir"/"bin

echo ""
echo "+ DOMINO Installation is FINISHED..."
echo ""
echo ""
echo ""
echo "+ Close the terminal and launch DOMINO using the Desktop link generated or using the command-line"

) 2>>DOMINO_Error.log