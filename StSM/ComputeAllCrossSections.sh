#This script is tested in CMSSW_5_0_0_pre5
#!/bin/bash
GlobalTag=MC_44_V1
if [ ! -d "$CMSSW_BASE/src/Configuration/GenProduction" ]; then
    mkdir $CMSSW_BASE/src/Configuration
    mkdir $CMSSW_BASE/src/Configuration/GenProduction
    mkdir $CMSSW_BASE/src/Configuration/GenProduction/python
fi
for ep in $( seq 20 10 20 )
do
    if [ -n "$2" ]&&[ "0.0$ep" != "$2" ]; then
	continue
    fi
    for mass in $( seq 200 25 1500 )
    do
	if [ -n "$1" ]&&[ "$mass" != "$1" ]; then
	    continue
	fi
	xsec_result=zprime_stu_m"$mass"_ep0"$ep"_cteq6l1_all.txt
	width_result=zprime_stu_m"$mass"_ep0"$ep"_cteq6l1_decaywidth.txt
	Written=FALSE
	CNAME=BatchMacro_"$mass"_"$ep"
	ShName=BatchMacro_M"$mass"_E0"$ep".sh
	card=StuZprimeToMuMu_M"$mass"_Epsilon0"$ep"_7TeV_pythia6_cff.py
	rootfile=/tmp/zhangjin/StuZprimeToMuMu_M"$mass"_Epsilon0"$ep"_7TeV_pythia6.root
	echo "#include \"TSystem.h\"
               void $CNAME () {
                 gSystem->Load(\"CouplingCons_C.so\");
                 Cal($mass,0.0$ep);
                 gApplication->Terminate(0);
               }">$CNAME.C
	root -l -b $CNAME.C
	rm $CNAME.C
	echo "#include \"TSystem.h\"
               void $CNAME () {
                 gSystem->Load(\"FWLite_GenShape_C.so\");
                 FWLite_GenShape(\"$rootfile\");
                 gApplication->Terminate(0);
               }">$CNAME.C
	echo "#include \"TSystem.h\"
               void final_$CNAME () {
                 gSystem->Load(\"WidthFit_C.so\");
                 WidthFit(\"${rootfile/.root/_GenShape.root}\");
                 gApplication->Terminate(0);
               }">final_$CNAME.C
	echo "#!/bin/bash
               cd `pwd`
               mv $card $CMSSW_BASE/src/Configuration/GenProduction/python
               cd $CMSSW_BASE/src/
               eval \`scramv1 runtime -sh\`
               sed -i \"s/crossSection.*cms.untracked.double(.*),/crossSection = cms.untracked.double(0),/g\" Configuration/GenProduction/python/$card
               cd Configuration/GenProduction
               scramv1 b
               cd $CMSSW_BASE/src/
               cmsDriver.py Configuration/GenProduction/python/$card -s GEN --fileout $rootfile --conditions $GlobalTag::All --datatier GEN --eventcontent RAWSIM -n 5000 --mc --no_exec
               echo 'process.load(\"FWCore.MessageLogger.MessageLogger_cfi\")'>>${card%.py}_py_GEN.py
               echo 'process.MessageLogger.cerr.FwkReport.reportEvery = 1000'>>${card%.py}_py_GEN.py
               CROSSSECTION=\`cmsRun ${card%.py}_py_GEN.py | grep \"All included\"\`
               rm ${card%.py}_py_GEN.py
               cd `pwd`
               if [ -n \"\$CROSSSECTION\" ]; then
                   Exponent=\`expr match \"\$CROSSSECTION\" \".*D\\([0-9+-]*\\)\"\`
                   CROSSSECTION=\`expr match \"\$CROSSSECTION\" \".*\\([0-9].[0-9]\\{3\\}\\)\"\`
                   echo \"$mass \$CROSSSECTION\"\"E\"\"\$Exponent\">$xsec_result
               fi
               root -l -b $CNAME.C
               root -l -b final_$CNAME.C>$CNAME.txt
               tail -n 1 $CNAME.txt>$width_result
               rm $CNAME.txt
               rm $CNAME.C
               rm final_$CNAME.C
               rm $rootfile
               rm ${rootfile/.root/_GenShape.root}
               rm $ShName
               ">$ShName
	chmod +x $ShName
        bsub -o /dev/null -e /dev/null -q 1nd $ShName
    done
done
