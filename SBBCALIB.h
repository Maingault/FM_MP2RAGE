



#ifndef __SBBCALIB_H
#define __SBBCALIB_H 1

// * -------------------------------------------------------------------------- *
// * Includes                                                                   *
// * -------------------------------------------------------------------------- *
#include "MrServers/MrImaging/libSBB/SeqBuildBlock.h"

#include "MrServers/MrMeasSrv/SeqIF/sde_allincludes.h"
#include "MrServers/MrMeasSrv/SeqIF/libMES/SEQData.h"

//new includes
#include "MrServers/MrImaging/libSeqSysProp/SysProperties.h"
#include "MrServers/MrProtSrv/MrProt/MrSliceGroup.h"
#include "MrServers/MrProtSrv/MrProt/MrSlice.h"
#include "MrServers/MrProtSrv/MrProt/MrCoilInfo.h"
#include "MrServers/MrProtSrv/MrProt/CoilSelect/MrCoilSelect.h"
#include "ProtBasic/Interfaces/MrProtData.h"
#include "ProtBasic/Interfaces/MrPat.h"
#include "ProtBasic/Interfaces/MrWipMemBlock.h"
#include "ProtBasic/Interfaces/MrProtWorkflow.h"
#include "MrServers/MrMeasSrv/SeqIF/libRT/sRF_PULSE.h"
#include "MrServers/MrMeasSrv/SeqIF/libRT/sFREQ_PHASE.h"
#include "MrServers/MrMeasSrv/SeqIF/libRT/sGRAD_PULSE.h"
//#include "MrServers/MrMeasSrv/SeqIF/libRT/sGRAD_BASE.h"
#include "MrServers/MrMeasSrv/SeqIF/libRT/sREADOUT.h"
//required for SliceInfo sending
#include "MrServers/MrMeasSrv/SeqIF/libRT/sSYNCDATA.h"
#include "ProtBasic/Interfaces/ExternalInterface.h"
#include "ProtBasic/Interfaces/MrTherapyExchangeData.h"


// * -------------------------------------------------------------------------- *
//  CALIBWC Class
// * -------------------------------------------------------------------------- *
class SBBCALIB: public SeqBuildBlock
{
public:
		
		//Constructor & Destructor : 
		SBBCALIB(SBBList* pSBBList);
		virtual ~SBBCALIB();
	
		virtual void getNbpoints(long Nbpoints);

		virtual void getRadial(bool m_bRadial);
	
		//  --------------------------------------------------------------------------
		//
		//  Name        :  SBBCALIB::prep
		//
		//  Description :
		/// \brief <b>     Preparation of the sequence during binary search and prior
		///                 to sequence execution  </b>
		//
		//  Return      :  NLS status
		//
		//  --------------------------------------------------------------------------
		
		virtual bool prep (MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo);
	
	
		
	
		//  --------------------------------------------------------------------------
		//
		//  Name        :  SBBCALIB::run
		//
 		//  Description :
		/// \brief <b>     Run the calibration of the Wave Caipi sequence </b>
		//  
		//  Return      :  NLS status
		//
		//  --------------------------------------------------------------------------
	
		virtual bool run (MrProt &rMrProt, SeqLim &rSeqLim, SeqExpo &rSeqExpo, sSLICE_POS* pSLC );



		//  --------------------------------------------------------------------------
		//
		//  Name        :  SBBCALIB::runKernel
		//
		//  Description :
		/// \brief <b>     Run the basic kernel of the Wave Caipi sequence </b>
		//  
		//  Return      :  NLS status
		//
		//  --------------------------------------------------------------------------
		bool runKernel  (MrProt &rMrProt,SeqLim &rSeqLim, SeqExpo &rSeqExpo, long lKernelMode,long partition ,long lSlice, long lLine);



		//  --------------------------------------------------------------------------
		// PUBLIC METHODS
		//  --------------------------------------------------------------------------
		long   getKernelRequestsPerMeasurement(void);

		long   getDurationMainEventBlock(void);
		
protected:
		

		//  ------------------------------------------------------------------
		///
		//  ------------------------------------------------------------------
		long       m_alTEFil[K_NO_TIME_ELEMENTS]               ;
		int32_t    m_lLinesToMeasure                           ;
		long       m_lRepetitionsToMeasure                     ;
		long       m_lPhasesToMeasure                          ;
		long       m_lSlicesToMeasure                          ;
		int32_t    m_lPartitionsToMeasure                      ;
		long       m_lScanTimeAllSats                          ;
		long       m_lDurationMainEventBlock                   ;
		long       m_lKernelRequestsPerMeasurement             ;
		long       m_lKernelCallsPerRelevantSignal             ;
		long       m_lInterDuration                            ;
		long       m_lTrigHaltDuration1                        ;
		long       m_lTrigHaltDuration2                        ;

		double     m_dRFSpoilIncrement                         ;
		double     m_dRFSpoilPhase                             ;
		long       m_lMySliSelRampTime                         ;
		double     m_dMinRiseTime                              ;
		double     m_dGradMaxAmpl                              ;


		// basic timing elements (Tot = total time)
		long         m_lExcTime                                ;                // begin sequence -> end RF-pulse            (Excitation Time)
		long         m_lEncTime                                ;                // end RF-pulse   -> begin ADC               (Encoding Time)
		long         m_lSampTime                               ;                // begin 1st ADC  -> end last ADC            (Sampling Time)
		long         m_lSpoilTime                              ;                // end ADC        -> end TR                  (Spoiling Time)                                                                                    (Spoiling/Unwinding)
		long         m_lTRMin                                  ;                // minimum TR
	    long         m_lStartSRFTime                           ;
        long         m_lSliSelRampTime                         ;

		// Slice position information (rotation matrices and shifts)
		sSLICE_POS m_asSLC[K_NO_SLI_MAX]                       ;


		// These variables are needed to remember the last excited slice to use the correct rf spoil phase
		double      m_dRFSpoilPhasePrevSlice          ;       // Remember RF spoil phase of previous slice
		double      m_dRFSpoilIncrementPrevSlice      ;       // Remember RF spoil phase increment
		double      m_dRFPrevSlicePosSag              ;  // These are used to check whether we're dealing
		double      m_dRFPrevSlicePosCor              ;  // with a different slice when RF spoiling is used
		double      m_dRFPrevSlicePosTra              ;
		double      m_dRFPrevSliceNormalSag           ;
		double      m_dRFPrevSliceNormalCor           ;
		double      m_dRFPrevSliceNormalTra           ;
       
		long m_lRCColumn;
		bool m_bRadial;
		long m_lSBBTime;

		//. ----------------------------------------------------------
		//. Instantiate RF Pulse objects
		//. ----------------------------------------------------------

		// every rf-pulse must (!) have an unique name
		//  (e.g. fl_Flash_ex = flash template excitation)"
		//  maximium of 12 chars
		sRF_PULSE_SINC           m_sSRF01                   ;
		// define an event that sets the transmitter phase
		sFREQ_PHASE              m_sSRF01zSet               ;
		// define an event that resets the transmitter phase
		sFREQ_PHASE              m_sSRF01zNeg              ;


		//. ----------------------------------------------------------
		//. Instantiate Gradient Pulse objects
		//. ----------------------------------------------------------

		// every gradient pulse event should have an unique name,
		// (e.g. the name of its structure)
		sGRAD_PULSE            m_sGSliSel                ;           // Slice-select pulse
		sGRAD_PULSE            m_sGSliSelReph            ;   // Slice-select rephaser pulse
	    // note that one can use arrays of gradient
		sGRAD_PULSE       m_sGSpoilSlice;
		sGRAD_PULSE       m_sGSpoilRead;

		//FMUTE
		sGRAD_PULSE m_sGradReadReph							;
		sGRAD_PULSE m_sGradPhaseReph						;
		sGRAD_PULSE m_sGradSliceReph						;
		sGRAD_PULSE m_sGradReadDephFM						;
		sGRAD_PULSE m_sGradReadRephFM						;


		// the readout gradient pulse event belongs to a different class:
		sGRAD_PULSE m_sGradRead								;	
		sGRAD_PULSE m_sGradReadDeph							;

		//. ----------------------------------------------------------
		//. Instantiate Readout objects
		//. ----------------------------------------------------------
		// every readout event must have an unique identifier
		sREADOUT               m_sADC[20]                   ;
		// define an event that sets the receiver phase
		sFREQ_PHASE            m_sADCzSet                  ;
		// define an event that undoes the receiver phase
		sFREQ_PHASE            m_sADCzNeg                  ;


		//. ----------------------------------------------------------
		//. Instantiate Sync objects: Osc bit and triggering
		//. ----------------------------------------------------------

		sSYNC_OSC              m_sOscBit;
		sSYNC_PHYSIO1_HALT     m_sTriggerBit1;
		sSYNC_PHYSIO2_HALT     m_sTriggerBit2;


	    double m_pauseBTrep;
		double m_AverageNum;

};
#endif