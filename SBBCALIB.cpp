#include "MrServers/MrMeasSrv/SeqIF/MDH/MdhProxy.h"
#include "MrServers/MrMeasSrv/SeqIF/Sequence/sequmsg.h"
#include "MrServers/MrMeasSrv/SeqIF/libRT/libRT.h"
#include "MrServers/MrMeasSrv/SeqIF/csequence.h"
#include "MrServers/MrMeasSrv/SeqFW/libGSL/libGSL.h"



#include "MrServers/MrProtSrv/MrProt/KSpace/MrKSpace.h"    // needed for accessing k-space parameters in protocol 
#include "MrServers/MrProtSrv/MrProtocol/UILink/StdProtRes/StdProtRes.h"
#include "MrServers/MrProtSrv/MrProt/CoilSelect/MrRxCoilSelect.h"
#include "MrServers/MrProtSrv/MrProt/CoilSelect/MrCoilSelect.h"
#include "MrServers/MrProtSrv/MrProt/MeasParameter/MrSysSpec.h" // for SysProperties
#include "MrServers/MrProtSrv/MrProt/MeasParameter/MrRXSpec.h"
#include "MrServers/MrProtSrv/MrProt/MrSliceGroup.h"
#include "MrServers/MrProtSrv/MrProt/MrCoilInfo.h"

#include "MrServers/MrImaging/ut/libsequt.h"                // for mSEQTest
#include "ProtBasic/Interfaces/MrRXSpec.h"
#include "ProtBasic/Interfaces/MrWipMemBlock.h"
#include "ProtBasic/Interfaces/MrInteractive.h"
#include "ProtBasic/Interfaces/MrInteractiveRealtime.h"
#include "ProtBasic/Interfaces/ExternalInterface.h"
#include "ProtBasic/Interfaces/MeasAqcDefs.h"
#include "ProtBasic/Interfaces/MrFastImaging.h"
#include <iostream>
#include "SBBCALIB.h"

#ifndef SEQ_NAMESPACE
#error SEQ_NAMESPACE not defined
#endif

double u_lDistShift ;
double u_lAverageNb ;
extern long l_offset;



#define OnErrorPrintAndReturn(S,P) if (!MrSucceeded(S)) \
{ MRTRACE("Error from %s \n",P); return(S);}

//+ [ contructor **************************************************************
//+ ***************************************************************************
// with SBB-list
SBBCALIB::SBBCALIB(SBBList* pSBBList)
: SeqBuildBlock (pSBBList),	
m_sSRF01            ("sSRF01"),
m_sSRF01zSet        ("sSRF01zSet"),
m_sSRF01zNeg        ("sSRF01zNeg")	,
m_sGSliSel                                    ("m_sGSliSel"),// Slice-select pulse
// m_sGReadDeph                                  ("m_sGReadDeph"),
m_sGSliSelReph                                ("m_sGSliSelReph"), // Slice-select rephaser pulse
m_sGradRead                               ("m_sGradReadoutref"),
m_sGradReadDeph                                ("m_sGradReadoutr"),
m_sADCzSet                                 ("sADC01zSet"),
m_sADCzNeg                                ("sADC01zNeg"),
m_lRCColumn (1),
m_pauseBTrep		(3000000)

//FMUTE
, m_sGradReadReph								("m_sGradReadReph")
, m_sGradPhaseReph								("m_sGradPhaseReph")
, m_sGradSliceReph								("m_sGradSliceReph")
, m_sGradReadDephFM								("m_sGradReadDephFM")
, m_sGradReadRephFM								("sGradReadRephFM")

{



};


//+ [ destructor **************************************************************
//+ ***************************************************************************
SBBCALIB::~SBBCALIB()
{

}

void SBBCALIB::getNbpoints(long Nbpoints)
{
	m_lRCColumn=Nbpoints;
}

void SBBCALIB::getRadial(bool bRadial)
{
	m_bRadial = bRadial;
}

bool SBBCALIB::prep(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo)
{	

	//. --------------------------------------------------------
	//. Declaration of local variables
	//. --------------------------------------------------------
	NLS_STATUS  lStatus = SEQU__NORMAL;

	m_dMinRiseTime = 1.2 * SysProperties::getGradMinRiseTime(rMrProt.gradSpec().mode());
	m_dGradMaxAmpl = 0.8 * SysProperties::getGradMaxAmpl(rMrProt.gradSpec().mode());

	u_lDistShift=12;
	m_AverageNum=1;

	//. ---------------------------------------------------------------------------
	//. Prepare the RF pulse objects
	//. ---------------------------------------------------------------------------

	// initialize variables for rf spoiler

	m_sSRF01.setFlipAngle            (90);
	m_sSRF01.setInitialPhase         (0);
	m_sSRF01.setThickness            (2);
	m_sSRF01.setTypeExcitation       ();
	m_sSRF01.setAsymmetry             (0.5);
	m_sSRF01.setDuration             (2560);
	// this it required by the unit test, so it knows the purpose the pulse
	m_sSRF01.setSamples              (128);
	m_sSRF01.setBandwidthTimeProduct (2.70);

	if (! m_sSRF01.prepSinc(rMrProt,rSeqExpo)) return (m_sSRF01.getNLSStatus());

	//. -------------------------------------------------------------------------------
	//. Prepare readout (ADC) objects
	//. -------------------------------------------------------------------------------
	double dReadAsymmetry = 0;// 100% asymmetric
	long lColumns = (long) ((rMrProt.kSpace().baseResolution()+ 2*l_offset) * (dReadAsymmetry + 0.5)  + 0.5);
	
	lColumns = m_lRCColumn;
	
	for(unsigned i=0; i<rMrProt.contrasts();i++)
		m_sADC[i].prep(lColumns, static_cast<long>(rMrProt.rxSpec().effDwellTime( rSeqLim.getReadoutOSFactor() )[0] ));

	//. -------------------------------------------------------------------------------
	//. Calculation of Gradient Amplitude
	//. -------------------------------------------------------------------------------
	double rampTime = 300;//TODO correct FOV depending on the rampsampling
	double GradientAmplitudeFOV= 1e9/(rMrProt.sliceSeries()[0].readoutFOV()*m_sADC[0].getDwellTime()*m_sSRF01.getLarmorConst());

	//. -------------------------------------------------------------------------------
	//. read phase gradient pulse																			a modifier en fonction de m_bRadial
	//. -------------------------------------------------------------------------------
	m_sGradRead.setMaxMagnitude(m_dGradMaxAmpl);
	m_sGradRead.setMinRiseTime (m_dMinRiseTime);
	m_sGradRead.setRampTimes(fSDSRoundUpGRT(rampTime));

	m_sGradRead.setDuration  (fSDSRoundUpGRT(m_sADC[0].getDuration())) ;  /*! EGA-03; EGA-01 !*/
	
	
	if(m_bRadial){
		m_sGradRead.setDuration  (fSDSRoundUpGRT(2*m_sGradRead.getDuration()-rampTime/2)) ;  /*! EGA-03; EGA-01 !*/
	}

	m_sGradRead.prepAmplitude(GradientAmplitudeFOV);
	
	if (!m_sGradRead.check()) return (m_sGradRead.getNLSStatus());

	//. -------------------------------------------------------------------------------
	//. read dephase gradient pulse																	FM
	//. -------------------------------------------------------------------------------
	if(m_bRadial){
		m_sGradReadDephFM.setMaxMagnitude(m_dGradMaxAmpl);
		m_sGradReadDephFM.setMinRiseTime (m_dMinRiseTime);
		m_sGradReadDephFM.setRampTimes(fSDSRoundUpGRT(rampTime));

		m_sGradReadDephFM.setDuration (fSDSRoundUpGRT(0.5*(m_sGradRead.getDuration() - 0.5*rampTime) ));
	
		std::cout << "durée de dephasage : " << m_sGradReadDephFM.getDuration() << std::endl;
		m_sGradReadDephFM.prepAmplitude(-GradientAmplitudeFOV);
	
		if (!m_sGradReadDephFM.check()) return (m_sGradReadDephFM.getNLSStatus());
	}

	//. -------------------------------------------------------------------------------
	//. read rephase gradient pulse																	FM
	//. -------------------------------------------------------------------------------
	if(m_bRadial){
		m_sGradReadRephFM.setMaxMagnitude(m_dGradMaxAmpl);
		m_sGradReadRephFM.setMinRiseTime (m_dMinRiseTime);
		m_sGradReadRephFM.setRampTimes(fSDSRoundUpGRT(rampTime));

	
		m_sGradReadRephFM.setDuration (fSDSRoundUpGRT(m_sGradRead.getDuration() - m_sGradReadDephFM.getDuration() ));
	
		std::cout << "durée de dephasage : " << m_sGradReadDephFM.getDuration() << std::endl;
		m_sGradReadRephFM.prepAmplitude(-GradientAmplitudeFOV);
	
		if (!m_sGradReadRephFM.check()) return (m_sGradReadRephFM.getNLSStatus());
	}


	//. -------------------------------------------------------------------------------
	//. read phase gradient pulse
	//. -------------------------------------------------------------------------------
	m_sGradReadDeph.setMaxMagnitude(m_dGradMaxAmpl);
	m_sGradReadDeph.setMinRiseTime (m_dMinRiseTime);
	m_sGradReadDeph.setRampTimes(fSDSRoundUpGRT(rampTime));
	m_sGradReadDeph.setDuration  (fSDSRoundUpGRT(m_sADC[0].getDuration())) ;             /*! EGA-03; EGA-01 !*/
	m_sGradReadDeph.prepAmplitude(-GradientAmplitudeFOV);

	if (!m_sGradReadDeph.check()) return (m_sGradReadDeph.getNLSStatus());

	//. ----------------------------------------------------------------------------
	//. Prepare slice selection
	//. ----------------------------------------------------------------------------	
	double dMomentumToBeRewinded = 0;
	m_lSliSelRampTime =  fSDSRoundUpGRT( std::max<long>(SysProperties::getCoilCtrlLead(), m_dMinRiseTime * m_sSRF01.getGSAmplitude()));
	dMomentumToBeRewinded =  m_sSRF01.getGSAmplitude() * ( m_sSRF01.getDuration() + m_lSliSelRampTime ) / 2.;
	if (!m_sGSliSel.prepAmplitude(m_lSliSelRampTime, fSDSRoundUpGRT(m_sSRF01.getDuration()+m_lSliSelRampTime), m_lSliSelRampTime, m_sSRF01.getGSAmplitude() ) ) /*! EGA-04; EGA-02 !*/
		return(m_sGSliSel.getNLSStatus());
	if (! m_sGSliSel.check() ) return (m_sGSliSel.getNLSStatus());
	m_sGSliSelReph.setMaxMagnitude(m_dGradMaxAmpl);
	m_sGSliSelReph.setMinRiseTime(m_dMinRiseTime);
	if (! m_sGSliSelReph.prepSymmetricTOTShortestTime(-dMomentumToBeRewinded))
		return m_sGSliSelReph.getNLSStatus();

	m_lStartSRFTime = m_lSliSelRampTime;
	m_lExcTime = m_lStartSRFTime + m_sSRF01.getDuration();
	m_lEncTime = m_sGSliSel.getRampDownTime() + m_sGSliSelReph.getTotalTime();


	//. -------------------------------------------------------------------------------  FMUTE
	//. read rephasing gradient pulse
	//. -------------------------------------------------------------------------------
	m_sGradReadReph.setMaxMagnitude(m_dGradMaxAmpl);
	m_sGradReadReph.setMinRiseTime (m_dMinRiseTime);
	m_sGradReadReph.setRampTimes(fSDSRoundUpGRT(rampTime));
	m_sGradReadReph.setDuration  (fSDSRoundUpGRT(m_sADC[0].getDuration())) ;             /*! EGA-03; EGA-01 !*/
	m_sGradReadReph.prepAmplitude(GradientAmplitudeFOV);

	m_sGradReadReph.setAmplitude(-m_sGradRead.getAmplitude());
	
	if (!m_sGradReadReph.check()) return (m_sGradReadReph.getNLSStatus());
	

	//. -------------------------------------------------------------------------------  FMUTE
	//. phase rephasing gradient pulse
	//. -------------------------------------------------------------------------------
	m_sGradPhaseReph.setMaxMagnitude(m_dGradMaxAmpl);
	m_sGradPhaseReph.setMinRiseTime (m_dMinRiseTime);
	m_sGradPhaseReph.setRampTimes(fSDSRoundUpGRT(rampTime));
	m_sGradPhaseReph.setDuration  (fSDSRoundUpGRT(m_sADC[0].getDuration())) ;             /*! EGA-03; EGA-01 !*/
	m_sGradPhaseReph.prepAmplitude(GradientAmplitudeFOV);
	
	if (!m_sGradPhaseReph.check()) return (m_sGradPhaseReph.getNLSStatus());
	

	//. -------------------------------------------------------------------------------  FMUTE
	//. slice rephasing gradient pulse
	//. -------------------------------------------------------------------------------
	m_sGradSliceReph.setMaxMagnitude(m_dGradMaxAmpl);
	m_sGradSliceReph.setMinRiseTime (m_dMinRiseTime);
	m_sGradSliceReph.setRampTimes(fSDSRoundUpGRT(rampTime));
	m_sGradSliceReph.setDuration  (fSDSRoundUpGRT(m_sADC[0].getDuration())) ;             /*! EGA-03; EGA-01 !*/
	m_sGradSliceReph.prepAmplitude(GradientAmplitudeFOV);
	
	if (!m_sGradSliceReph.check()) return (m_sGradSliceReph.getNLSStatus());

	if(m_bRadial)
	{
		lColumns = 3*lColumns;
	}
	for(unsigned i=0; i<rMrProt.contrasts();i++)
		m_sADC[i].prep(lColumns, static_cast<long>(rMrProt.rxSpec().effDwellTime( rSeqLim.getReadoutOSFactor() )[0] ));

	//. ----------------------------------------------------------------------------
	//. Calculate TEFill-times and check, whether timing can be realized
	//. ----------------------------------------------------------------------------
	double MinDurationBetweenADCandRF=1.2 * SysProperties::getMinDurationBetweenRFPulseAndReadout();
	double MinDurationBetweenRFandADC = 40;// fm_ute WIP SIEMENS

	m_lSBBTime = m_sGSliSel.getTotalTime() + m_sGSliSelReph.getTotalTime() + m_sGradRead.getTotalTime() + m_pauseBTrep;
	if(m_bRadial){
		m_lSBBTime += m_sGradReadRephFM.getTotalTime() + m_sGradReadDephFM.getTotalTime();
	}


	m_lSBBTime = 6 * m_lSBBTime;

	m_lSBBDurationPerRequest_us = m_lSBBTime;

	std::cout << "m_lSBBTime ===== "<< m_lSBBTime << std::endl;
	
	m_RFInfoPerRequest = m_sSRF01.getRFInfo()*6;

	for(unsigned i=1; i<rMrProt.contrasts();i++){

		m_alTEFil[ 0] = rMrProt.te()[0]  - (m_sSRF01.getDuration()/2 + MinDurationBetweenRFandADC);

		if ( m_alTEFil[ 0] < 0 )
			if (rSeqLim.isContextPrepForMrProtUpdate())
				rMrProt.te()[0] += fSDSRoundUpGRT( -m_alTEFil[ 0]);
			else
				return SBB_NEGATIV_TEFILL ;
		m_alTEFil[ 0] = fSDSRoundUpGRT(m_alTEFil[ 0]);


		m_alTEFil[i] = rMrProt.te()[i] - (rMrProt.te()[i-1] + 2 * m_sGradRead.getTotalTime());

		if ( m_alTEFil[ i] < 0 )
			if (rSeqLim.isContextPrepForMrProtUpdate())
				rMrProt.te()[i] += fSDSRoundUpGRT( -m_alTEFil[i]);
			else
				return SBB_NEGATIV_TEFILL ;
		m_alTEFil[ i] = fSDSRoundUpGRT(m_alTEFil[ i]);

	}

	lStatus = fSUPrepSlicePosArray (rMrProt, rSeqLim, m_asSLC);
	setPrepared();
}

bool SBBCALIB::run(MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo,  sSLICE_POS* pSLC )
{
	//. --------------------------------------------------------
	//. Local variables
	//. --------------------------------------------------------
	static const char *ptModule = {"fSEQRun"};
	NLS_STATUS lStatus          = SEQU__NORMAL;

	long       lRepetition       ;
	long       lLine             ;
	long       lPhase            ;
	long       lChronologicSlice ;

	//	mPrintTrace1 (DEBUG_RUN, DEBUG_CALL, "() <%s> started", rSeqLim.getLinkedSeqFilename() ) ;

	// Initialization of the unit test function
	// mSEQTest(rMrProt,rSeqLim,rSeqExpo,RTEB_ORIGIN_fSEQRunStart,0,0,0,0,0); /*! EGA-All !*/

	// Set the line independent MDH parameters
	m_sADC[0].getMDH().setKSpaceCentreLineNo      ( rMrProt.PAT().getlRefLinesPE()/2);
	if (rMrProt.kSpace().getucDimension() == SEQ::DIM_3)
	{
		m_sADC[0].getMDH().setKSpaceCentrePartitionNo (rMrProt.PAT().getlRefLinesPE()/2);
	}
	else
	{
		m_sADC[0].getMDH().setKSpaceCentrePartitionNo (0);
	}

	//. ---------------------------------------------------------------------------
	//. Measurement Repetitions
	//. ---------------------------------------------------------------------------
	for (lChronologicSlice = 0; lChronologicSlice < 6; lChronologicSlice++ )
	{
		for (lRepetition = 0; lRepetition < m_AverageNum  ; lRepetition++) // <= because measurements = repetitions + 1
		{	
			//mSEQTest (rMrProt, rSeqLim, rSeqExpo, RTEB_ClockInitTR, 0, 0, lChronologicSlice, 0, 0);

			//. ---------------------------------------------------------------------------
			//. Kernel Call
			//. ---------------------------------------------------------------------------
			lStatus = runKernel( rMrProt, rSeqLim, rSeqExpo, KERNEL_IMAGE, lChronologicSlice,0, lRepetition);
			OnErrorPrintAndReturn(lStatus,"runKernel");

		} //END of slice loop

	} //END of repetition loop

	// Tell sequence unit test that sequence run is finished.
	// mSEQTest(rMrProt,rSeqLim,rSeqExpo,RTEB_ORIGIN_fSEQRunFinish,0,0,0,0,0); /*! EGA-All !*/ 
	//	mPrintTrace1 (DEBUG_RUN, DEBUG_CALL | DEBUG_RETURN, "() <%s> finished", rSeqLim.getLinkedSeqFilename() ) ;

	return(lStatus);
}


bool SBBCALIB::runKernel  (MrProt &rMrProt,SeqLim &rSeqLim, SeqExpo &rSeqExpo, long lKernelMode, long lChronologicSlice,long partition, long lLine)
{

	//. --------------------------------------------------------------------------
	//. Local variables
	//. --------------------------------------------------------------------------
	static const char *ptModule         = {"runKernel"} ; // point to name of current module
	NLS_STATUS         lStatus          = SEQU__NORMAL ;      // a status variable
	unsigned long      ulTestIdent      = 0 ;                 // tell unit test whether we're running or checking the kernel
	long               lT               = 0 ;                 // used as clock time in the main event block

	if (lKernelMode == KERNEL_CHECK) ulTestIdent = RTEB_ORIGIN_fSEQCheck;
	else                             ulTestIdent = RTEB_ORIGIN_fSEQRunKernel;


	//. --------------------------------------------------------------------------
	//. Set the  MDH parameters
	//. --------------------------------------------------------------------------
	for(unsigned i=0; i< rMrProt.contrasts() ; i++){
		m_sADC[i].getMDH().setCslc (0); // what's the current slice
		m_sADC[i].getMDH().setCeco (i); // what's the current slice
		m_sADC[i].getMDH().setClin (lLine );             // what's the current line?
		m_sADC[i].getMDH().setCpar (lChronologicSlice);        // what's the current partition?
		
		m_sADC[i].getMDH().setEvalInfoMask (m_sADC[i].getMDH().getEvalInfoMask() | MDH_ONLINE) ;
		m_sADC[i].getMDH().setEvalInfoMask (m_sADC[i].getMDH().getEvalInfoMask() | MDH_PATREFSCAN) ;
	}


	switch(lChronologicSlice)
	{
	case 0:
		m_asSLC[0].setSliceShift(-u_lDistShift);
		break;
	case 1:
		m_asSLC[0].setSliceShift(u_lDistShift);
		break;
	case 2:
		m_asSLC[0].setSliceShift(-u_lDistShift);
		break;
	case 3:
		m_asSLC[0].setSliceShift(u_lDistShift);
		break;
	case 4:
		m_asSLC[0].setSliceShift(-u_lDistShift);
		break;
	case 5:
		m_asSLC[0].setSliceShift(u_lDistShift);
		break;
	}	

	//. --------------------------------------------------------------------------
	//. Set the frequency/phase properties of the RF pulses
	//. --------------------------------------------------------------------------
	m_sSRF01zSet.prepSet (m_asSLC[0], m_sSRF01) ;    /*! EGA-05 !*/
	m_sSRF01zNeg.prepNeg (m_asSLC[0], m_sSRF01) ;    /*! EGA-05 !*/


	//. --------------------------------------------------------------------------
	//. Set the frequency/phase properties of the adc OFFcenter is dealt in gadgetron
	//. --------------------------------------------------------------------------
	m_sADCzSet.set("sADC01zSet", 0, m_sSRF01.getInitialPhase());
	m_sADCzNeg.set("sADC01zNeg", 0, -m_sSRF01.getInitialPhase());


	//. ---------------------------------------------------------------------------
	//. Begin of event block
	//. ---------------------------------------------------------------------------
	fRTEBInit(m_asSLC[0].getROT_MATRIX());    /* EGA-07 */


	//- ***************************************************** S E Q U E N C E   T I M I N G *************************
	//- *           Start Time    |    NCO    |  SRF  |  ADC  |             Gradient Events             | Sync
	//- *             (usec)      |   Event   | Event | Event |     phase    |   read     |    slice    | Event
	//- *fRTEI(                   ,           ,       ,       ,              ,            ,             ,         );
	//- *************************************************************************************************************	
	lT = 0  ;
	fRTEI(m_lStartSRFTime                      ,&m_sSRF01zSet, &m_sSRF01,      0,             0,            0,            0,        0);

	fRTEI(m_lStartSRFTime+m_sSRF01.getDuration() ,&m_sSRF01zNeg,       0,      0,             0,            0,            0,        0);

	if ((lChronologicSlice == 0) || (lChronologicSlice == 1)) {
		fRTEI(                       lT,          0,       0,      0,         0    ,            &m_sGSliSel,    0,        0);	
		fRTEI(lT + m_sGSliSel.getTotalTime(),            0,      0,      0,          0  ,          &m_sGSliSelReph,0,        0);	

		lT += m_lEncTime + m_lExcTime;

		if(!m_bRadial){
			fRTEI(lT + fSDSRoundUpGRT((l_offset*m_sADC[0].getDwellTime())/(1000))   ,            0,        0,         0, 0  , &m_sGradRead  ,0   ,    0);
		}
		else{
			fRTEI(lT + fSDSRoundUpGRT((l_offset*m_sADC[0].getDwellTime())/(1000))   ,            0,        0,         0, 0  , &m_sGradReadDephFM  ,0   ,    0);
			fRTEI(lT + fSDSRoundUpGRT((l_offset*m_sADC[0].getDwellTime())/(1000)+ m_sGradReadDephFM.getTotalTime())   ,            0,        0,         0, 0  , &m_sGradRead  ,0   ,    0);
			fRTEI(lT + fSDSRoundUpGRT((l_offset*m_sADC[0].getDwellTime())/(1000)+ m_sGradReadDephFM.getTotalTime()) + m_sGradRead.getTotalTime()   ,            0,        0,         0, 0  , &m_sGradReadRephFM  ,0   ,    0);
		}

		fRTEI(lT                              ,&m_sADCzSet  ,        0, &m_sADC[0]  ,             0,            0,             0,    0);
		fRTEI(lT + m_sADC[0].getRoundedDuration(),&m_sADCzNeg  ,        0,         0,             0,            0,             0,    0);

		lT +=   m_sGradRead.getTotalTime() ;
		if(m_bRadial){
			lT += m_sGradReadDephFM.getTotalTime() + m_sGradReadRephFM.getTotalTime();
			//lT += 2*m_sGradReadDephFM.getTotalTime();
		}
	
		
		
		for(unsigned i=1; i< rMrProt.contrasts() ; i++){
			fRTEI(lT                          ,            0,        0,         0, 0  , &m_sGradReadDeph  ,0   ,    0);
			lT += m_sGradReadDeph.getTotalTime() + m_alTEFil[i];		

			fRTEI(lT                              ,            0,        0,         0,   0, &m_sGradRead  ,0   ,    0);

			fRTEI(lT                              ,&m_sADCzSet  ,        0, &m_sADC[i]  ,             0,            0,             0,    0);
			fRTEI(lT + m_sADC[i].getRoundedDuration(),&m_sADCzNeg  ,        0,         0,             0,            0,             0,    0);
			lT +=   m_sGradRead.getTotalTime() ;
			fRTEI(lT + m_sADC[i].getRoundedDuration(),&m_sADCzNeg  ,        0,         0,             0,            0,             0,    0);
			
		
		}
		
	}



	if ((lChronologicSlice == 2) || (lChronologicSlice ==3)) {
		fRTEI(                       lT ,          0,       0,      0,             &m_sGSliSel,            0,    0,        0);	
		fRTEI(lT + m_sGSliSel.getTotalTime(),            0,      0,      0,             &m_sGSliSelReph,            0,0,        0);	

		lT += m_lEncTime + m_lExcTime;

		if(!m_bRadial){	
			fRTEI(lT  + fSDSRoundUpGRT((l_offset*m_sADC[0].getDwellTime())/(1000))   ,            0,        0,         0, &m_sGradRead  , 0  ,0   ,    0);
		}
		else{
			fRTEI(lT + fSDSRoundUpGRT((l_offset*m_sADC[0].getDwellTime())/(1000))   ,            0,        0,         0, &m_sGradReadDephFM  , 0  ,0   ,    0);
			fRTEI(lT + fSDSRoundUpGRT((l_offset*m_sADC[0].getDwellTime())/(1000)+ m_sGradReadDephFM.getTotalTime())   ,            0,        0,         0, &m_sGradRead  , 0  ,0   ,    0);
			fRTEI(lT + fSDSRoundUpGRT((l_offset*m_sADC[0].getDwellTime())/(1000)+ m_sGradReadDephFM.getTotalTime()+m_sGradRead.getTotalTime())   ,            0,        0,         0, &m_sGradReadRephFM  , 0  ,0   ,    0);
		}

		
		fRTEI(lT                              ,&m_sADCzSet  ,        0, &m_sADC[0]  ,             0,            0,             0,    0);
		fRTEI(lT + m_sADC[0].getRoundedDuration(),&m_sADCzNeg  ,        0,         0,             0,            0,             0,    0);

		lT +=   m_sGradRead.getTotalTime() ;
		
		if(m_bRadial){
			lT += m_sGradReadDephFM.getTotalTime() + m_sGradReadRephFM.getTotalTime();
			//lT += 2*m_sGradReadDephFM.getTotalTime();
		}
	
		for(unsigned i=1; i< rMrProt.contrasts() ; i++){
			fRTEI(lT                          ,            0,        0,         0, &m_sGradReadDeph  , 0  ,0   ,    0);
			lT += m_sGradReadDeph.getTotalTime() + m_alTEFil[i];		

			fRTEI(lT                              ,            0,        0,         0, &m_sGradRead  , 0  ,0   ,    0);

			fRTEI(lT                              ,&m_sADCzSet  ,        0, &m_sADC[i]  ,             0,            0,             0,    0);
			fRTEI(lT + m_sADC[i].getRoundedDuration(),&m_sADCzNeg  ,        0,         0,             0,            0,             0,    0);
			lT +=   m_sGradRead.getTotalTime() ;
		
		}
	}


	if ((lChronologicSlice == 4) || (lChronologicSlice == 5)) {

		fRTEI(                       lT ,          0,       0,      0,             0,            0,    &m_sGSliSel,        0);	
		fRTEI(lT + m_sGSliSel.getTotalTime(),            0,      0,      0,             0,            0,&m_sGSliSelReph,        0);	

		lT += m_lEncTime + m_lExcTime;

		if(!m_bRadial){
			fRTEI(lT   + fSDSRoundUpGRT((l_offset*m_sADC[0].getDwellTime())/(1000))  ,            0,        0,         0,  0 ,  0 ,&m_sGradRead   ,    0);
		}
		else{
			fRTEI(lT + fSDSRoundUpGRT((l_offset*m_sADC[0].getDwellTime())/(1000))   ,            0,        0,         0, 0  , 0  ,&m_sGradReadDephFM  ,    0);
			fRTEI(lT + fSDSRoundUpGRT((l_offset*m_sADC[0].getDwellTime())/(1000)+ m_sGradReadDephFM.getTotalTime())   ,            0,        0,         0, 0  , 0  ,&m_sGradRead   ,    0);
			fRTEI(lT + fSDSRoundUpGRT((l_offset*m_sADC[0].getDwellTime())/(1000)+ m_sGradRead.getTotalTime()+ m_sGradReadDephFM.getTotalTime())   ,            0,        0,         0, 0  , 0  ,&m_sGradReadRephFM   ,    0);
		}


		fRTEI(lT                              ,&m_sADCzSet  ,        0, &m_sADC[0]  ,             0,            0,             0,    0);
		fRTEI(lT + m_sADC[0].getRoundedDuration(),&m_sADCzNeg  ,        0,         0,             0,            0,             0,    0);

		lT +=   m_sGradRead.getTotalTime() ;
		if(m_bRadial){
			lT += m_sGradReadDephFM.getTotalTime() + m_sGradReadRephFM.getTotalTime();
			//lT += 2*m_sGradReadDephFM.getTotalTime();
		}
	
		

		for(unsigned i=1; i< rMrProt.contrasts() ; i++){
			fRTEI(lT                          ,            0,        0,         0,  0 , 0  ,&m_sGradReadDeph   ,    0);
			lT += m_sGradReadDeph.getTotalTime() + m_alTEFil[i];		

			fRTEI(lT                              ,            0,        0,         0, 0  , 0  ,&m_sGradRead   ,    0);

			fRTEI(lT                              ,&m_sADCzSet  ,        0, &m_sADC[i]  ,             0,            0,             0,    0);
			fRTEI(lT + m_sADC[i].getRoundedDuration(),&m_sADCzNeg  ,        0,         0,             0,            0,             0,    0);
			lT +=   m_sGradRead.getTotalTime() ;
			
		}
	}

	fRTEI(lT + m_pauseBTrep ,   0,   0,         0,             0,            0,             0,    0);
	

	//mSEQTest (rMrProt, rSeqLim, rSeqExpo, ulTestIdent    , 10, lLine, m_asSLC[lChronologicSlice].getSliceIndex(), 0, 0) ;
	//mSEQTest (rMrProt, rSeqLim, rSeqExpo, RTEB_ClockCheck, 10, lLine, m_asSLC[lChronologicSlice].getSliceIndex(), 0, 0) ;
	OnErrorPrintAndReturn(lStatus = fRTEBFinish(),"fRTEBFinish [*0100*]");
	//. ---------------------------------------------------------------------------
	//. End of event block
	//. ---------------------------------------------------------------------------
	return(lStatus);

}

