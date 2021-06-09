#!/bin/bash

#Enter path to directory with the 'scripts' folder in it. Ex: FOLDER = '/home/dhamelin/CoVescape/'
FOLDER='/Users/davidhamelin/Documents/Graduate_schoold_stuff/COVID_MUT_PIPELINE_TEST/MHC_Class_I/COMBINE_ALL_CODE_FOR_GITHUB_CELLSYSTEMS/TEST4_AddMutationalPatterns'

####################################### DO NOT CHANGE CODE BELOW THIS LINE

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -h|--help)
    echo "--script RunNetMHCpan.sh"
    echo "inputs:"
    echo "  -r | --RawMutations: File of raw mutations in .txt format."
    echo "  -h | --HLAlist: list of class I HLA molecules to process. Should provide them in the following format: HLA-A01:01,HLA-A02:01,..."
    echo "  -o | --OutputFilePATH: path/to/output folder in which to store the program outputs"
    echo "  -p | --netMHCpanPATH: path/to/netMHCpan shell script"
    echo ""
    echo "--script PostnetMHC.sh"
    echo "inputs:"
    echo "  -n | --NetMHCOutput: netMHCpan output."
    echo "  -m | --MutationList: Formatted and translated mutation list"
    echo ""
    echo "--script GetCores.sh"
    echo "inputs:"
    echo "  -n | --NetMHCOutput: netMHCOutput"
    echo "  -t | --tempFrame: HLA-specific processed netMHCpan output. Includes all mutations which lead to a strong gain or loss of binding (%Rank>2 <--> %Rank<0.5), as well as all related information (protein, mutated/reference sequences, peptide length, %Rank FC)"
    echo ""
    echo "--script PetersFigure.sh"
    echo "inputs:"
    echo "  -c | --CoreOutput: HLA-specific processed netMHCpan output with information regarding residues predicted to directly interact with binding groove."
    echo ""
    echo "--script MutationalPatterns.sh"
    echo "  -m | --MutationList: Formatted and translated mutation list"
    

    

    #EXTENSION="$2"
    shift # past argument
    #shift # past value
    ;;
    -s|--script)
    Script="$2"
    shift # past argument
    shift # past value
    ;;

    #########################################RunnetMHCpan
    -r|--RawMutations)
    mutFile1="$2"
    shift # past argument
    shift # past value
    ;;
    -h|--HLAlist)
    HLAs="$2"
    shift # past argument
    shift # past value
    ;;

    -o|--OutputFilePATH)
    outputFilePATH="$2"
    shift # past argument
    shift # past value
    ;;

    -p|--netMHCpanPATH)
    netMHCpanPATH="$2"
    shift # past argument
    shift # past value
    ;;


    ########################################PostNetMHCpan
    -n|--NetMHCOutput)
    netMHCout="$2"
    shift # past argument
    shift # past value
    ;;

    -m|--MutationList)
    mutations="$2"
    shift # past argument
    shift # past value
    ;;


    ##############################################GetCores
    -t|--tempFrame)
    tempframe="$2"
    shift # past argument
    shift # past value
    ;;

    ##########################################PetersFigure
    -c|--CoreOutput)
    coreOut="$2"
    shift # past argument
    shift # past value
    ;;







    --default)
    DEFAULT=YES
    shift # past argument
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

#echo "Script  = ${Script}"
#echo "First word     = ${WordOne}"
#echo "Second Word    = ${WordTwo}"
#echo "DEFAULT         = ${DEFAULT}"
#echo "Number files in SEARCH PATH with EXTENSION:" $(ls -1 "${SEARCHPATH}"/*."${EXTENSION}" | wc -l)
#if [[ -n $1 ]]; then
#    echo "Last line of file specified as non-opt/last argument:"
#    tail -1 "$1"
#fi

if [ ! -z "${Script}"  ] ; then




    if [ $Script == "RunNetMHCpan.sh" ]; then
        echo "Running RunNetMHCpan.sh"
        if [ ! -z "${mutFile1}"  ] && [ ! -z "${HLAs}" ] && [ ! -z "${outputFilePATH}" ] && [ ! -z "${netMHCpanPATH}" ]; then
            echo "We're here"
            echo "${mutFile1}"

            $FOLDER/scripts/RunNetMHC/$Script $mutFile1 $HLAs $outputFilePATH $netMHCpanPATH $FOLDER
        else
            echo "one or more arguments are missing"

        fi



    elif [ $Script == "PostnetMHC.sh" ]; then
        echo "Running PostnetMHC.sh"
        if [ ! -z "${netMHCout}"  ] && [ ! -z "${mutations}" ]; then
            echo "We're here"
            $FOLDER/scripts/PostNetMHCPan/$Script $netMHCout $mutations $FOLDER

        else
            echo "one or more arguments are missing"

        fi

    


    elif [ $Script == "GetCores.sh" ]; then
        echo "Running GetCores.sh"
        if [ ! -z "${netMHCout}" ] && [ ! -z "${tempframe}" ]; then
            echo "We're here"
            $FOLDER/scripts/GetCore/$Script $netMHCout $tempframe $FOLDER

        else
            echo "one or more arguments are missing"

        fi

    


    elif [ $Script == "PetersFigure.sh" ]; then
        echo "Running PetersFigure.sh"
        if [ ! -z "${coreOut}" ]; then
            echo "We're here"
            $FOLDER/scripts/PetersFigure/$Script $coreOut $FOLDER

        else
            echo "one or more arguments are missing"

        fi



    elif [ $Script == "MutationalPatterns.sh" ]; then
        echo "Running MutationalPatterns.sh"
        if [ ! -z "${mutations}" ]; then
            echo "We're here"
            $FOLDER/scripts/MutationalPatterns/$Script $mutations $FOLDER 

        else
            echo "one or more arguments are missing"

        fi



    else echo "option does not correspond to any script."

    fi

else
    echo ""
    echo ""
    echo "No scripts were provided"

fi










