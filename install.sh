#! /bin/sh
(
echo "\n\n##########################################"
echo "# Executing DOMINO installation script:"
echo "##########################################\n\n"
echo "This script DOES NOT install the DOMINO GUI interface, only DOMINO perl scripts, binaries and"
echo "modules would be installed. Please use the installers already provided in the DOMINO_installers"
echo "folders for whole installation: GUI, DOMINO perl scripts, binaries, etc...\n\n"

echo "**** REQUIREMENTS ****"
echo "Please make sure you accomplish all the requirements:"
echo " - Perl Programming Language version (> v5.2)"
echo " - C++ compiler"
echo "	- Mac OS users may need to install Xcode package via App Store"
echo "	- For Linux users install build essential package o g++ compiler"
echo "		Ubuntu/Linux Mint:"
echo "			sudo apt-get install build-essential"
echo "			or "
echo "			sudo apt-get install g++"
echo ""
echo " - zlib compression library (samtools requirement) <http://www.zlib.net/>"
echo "	- Mac OS: already included and installed."
echo "	- For Linux users install zlib1g-dev package"
echo "		Ubuntu/Linux Mint:"
echo "			sudo apt-get install zlib1g-dev"
echo ""
echo ""
installer=$0
current_dir=$(dirname $installer)
install_dir=`pwd`

if [ $current_dir = "." ] 
then
	current_dir=`pwd`
fi

echo "Do you want to proceed with the DOMINO perl script and module installation?\n"
echo "Type 'yes' and press [ENTER]"
echo "Type 'no' and press [ENTER]\n"
read input

if [ $input = "no" ]
then 
	echo "\n\nThank you very much, please fulfil installation requirements "
	echo "or check DOMINO_installer folder and double click in the desired installer...\n"
	exit
fi

echo ""
echo "Please enter wether Linux or Mac, for installation system...\n" 
read OpSys
echo ""
echo "\n\n+ Checking parameters:"
## Checking option
if [ $OpSys = Linux ]
then
	echo "\t+ Launching linux installer...\n"
	binaries_tar=$current_dir"/files/"third_party_binaries-src.tar-Linux.gz
	file_name=third_party_binaries-src.tar-Linux.gz
elif [ $OpSys = Mac ]
then
	echo "\t+ Launching Mac installer...\n"
	binaries_tar=$current_dir"/files/"third_party_binaries-src.tar-Mac.gz
	file_name=third_party_binaries-src.tar-Mac.gz
else
	echo "\n+ Neither Linux or Mac string provided!!! \n\nAborting the installation...\n"
	exit
fi

echo "\n\n+ Launching the installer...\n"
echo "+ Working dir: "$current_dir
echo "+ Installation dir:"$install_dir
echo ""
scripts_folder=$install_dir"/bin"
mkdir $scripts_folder
echo "+ Generate a folder for DOMINO bin: "$scripts_folder
tmp_dir_name="tmp_dir"
tmp_dir_path=$install_dir"/"$tmp_dir_name
echo ""
echo "+ Generate an intermediate folder: "$tmp_dir_path
mkdir $tmp_dir_path
echo ""
echo "+ Copy, extract and compile necessary files..."
echo ""
cp $binaries_tar $tmp_dir_path
cd $tmp_dir_path
tar -zxvf "./"$file_name
echo ""
cd binaries
ls | while read folder_tar_gz
do
	echo "- Extracting and compiling: "$folder_tar_gz
	cd $folder_tar_gz
	ls | while read tar_gz
	do
		tar -zxvf $tar_gz
		rm $tar_gz

		if [ "$folder_tar_gz" = "SPADES" ];
		then
			if [ $OpSys = "Linux" ]
			then
				mkdir -p $scripts_folder"/SPAdes-3.8.1-Linux"
				ls "SPAdes-3.8.1-Linux/*" | while read files
				do
					cp -r "SPAdes-3.8.1-Linux/$files" $scripts_folder"/SPAdes-3.8.1-Linux"
				done
			fi
		fi

		if [ "$folder_tar_gz" = "CAP3" ];
		then
			mkdir -p $scripts_folder"/"cap3"/"bin
			if [ $OpSys = "Linux" ]
			then
				mv CAP3"/"cap3 $scripts_folder"/"cap3"/"bin
				mv CAP3"/"README $scripts_folder"/"cap3
				mv CAP3"/"doc $scripts_folder"/"cap3
			
			elif [ $OpSys = "Mac" ]
			then
				mv cap3.macosx.intel64"/"cap3 $scripts_folder"/"cap3"/"bin
				mv cap3.macosx.intel64"/"README $scripts_folder"/"cap3
				mv cap3.macosx.intel64"/"doc $scripts_folder"/"cap3
			fi
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
				if [ $OpSys = "Linux" ]
			then
				mira_folder=mira_4.0.2_linux-gnu_x86_64_static
				mv $mira_folder"/"bin"/"mira $scripts_folder"/"mira_v4.0"/"bin
				mv $mira_folder"/"LICENCE $scripts_folder"/"mira_v4.0
				mv $mira_folder"/"README $scripts_folder"/"mira_v4.0
			
			elif [ $OpSys = "Mac" ]
			then
				mira_folder=mira_4.0.2_darwin13.1.0_x86_64_static
				mv $mira_folder"/"bin"/"mira $scripts_folder"/"mira_v4.0"/"bin
				mv $mira_folder"/"LICENCE $scripts_folder"/"mira_v4.0
				mv $mira_folder"/"README $scripts_folder"/"mira_v4.0
				mv $mira_folder"/"lib $scripts_folder"/"mira_v4.0
			fi
		fi			

		if [ "$folder_tar_gz" = "NGSQC" ]
		then
			mkdir -p $scripts_folder"/"NGSQCToolkit_v2.3.1
			mv NGSQCToolkit_v2.3.3"/"* $scripts_folder"/"NGSQCToolkit_v2.3.1
		fi

		if [ "$folder_tar_gz" = "NCBI" ]
		then
			mkdir -p $scripts_folder"/"NCBI_BLAST_v2.2.28
			if [ $OpSys = "Linux" ]
			then
				cd ncbi-blast-2.2.28+-src"/"c++
				sh ./configure > DOMINO_Error.log
	 			make > DOMINO_Error.log
				mv *"/"*bin"/"blastn $scripts_folder"/"NCBI_BLAST_v2.2.28
				mv *"/"*bin"/"makeblastdb $scripts_folder"/"NCBI_BLAST_v2.2.28
			
			elif [ $OpSys = "Mac" ]
			then
				mv * $scripts_folder"/"NCBI_BLAST_v2.2.28
			fi
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
cd $current_dir
echo ""
echo "+ All necessary binaries have been extracted, compiled and move to its place..."
echo ""
echo "+ Compile some perl modules..."
cp $current_dir"/files/tarGZ_modules_Perl_NGSQCToolkit.tar.gz" $tmp_dir_path
cd $tmp_dir_path
tar -zxvf tarGZ_modules_Perl_NGSQCToolkit.tar.gz
cd $current_dir
perl $current_dir"/files/unzip_install_NGSQCtoolkit_modules.pl" $tmp_dir_path"/"tarGZ_modules_Perl_NGSQCToolkit $install_dir"/bin/NGSQCToolkit_v2.3.1"
cp $current_dir"/files/tarGZ_modules_Perl.tar.gz" $tmp_dir_path
cd $tmp_dir_path
tar -zxvf tarGZ_modules_Perl.tar.gz
cd $current_dir
perl $current_dir"/files/unzip_install_perl_modules.pl" $tmp_dir_path"/"tarGZ_modules_Perl $install_dir"/bin/"
cp -r $install_dir"/bin/NGSQCToolkit_v2.3.1/lib/Parallel" $install_dir"/bin/lib"

cp $current_dir"/files/db_default.tar.gz" $install_dir"/bin/"
cd $install_dir"/bin/"
tar -zxvf "db_default.tar.gz"
rm -r "db_default.tar.gz"
cd $install_dir
echo "+ Removing temporary files...\n"
rm -rf $tmp_dir_path

## Copy DOMINO perl code
ls "$current_dir/src/perl/" | while read files
do
	cp "$current_dir/src/perl/$files" $install_dir"/bin"
done

mv "$current_dir/src/perl/DOMINO.pm" $install_dir"/bin/lib"

) 2>>DOMINO_error.log
