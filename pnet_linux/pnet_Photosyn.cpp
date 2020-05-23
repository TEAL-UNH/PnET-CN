#include "pnet_model.h"

void pnet_model::Photosyn(site_struct* site, veg_struct* veg, clim_struct* clim, int rstep, share_struct* share)
{
	//
	// Photosynthesis routine for the PnET ecosystem model.
	//


	int  ix;
	double i;
	double PsnTMax, DTemp, ten9, Ca, CiCaRatio, Ci350, CiElev, Arel350, ArelElev, Delgs, gsSlope, gsInt, Amax, BaseFolResp; 
	double GrossAmax, CanopyNetPsnO3, O3Prof, Layer, SLWLayer, Il, LightEff, LayerGrossPsnRate, LayerGrossPsn;
	double LayerResp, LayerNetPsn, netPsnumol, RelLayer, RelO3, LayerG, LayerDO3, LayerNetPsnO3, CanopyNetPsnPot;

	PsnTMax = veg->PsnTOpt + (veg->PsnTOpt - veg->PsnTMin);
	DTemp = ((PsnTMax - share->Tday) * (share->Tday - veg->PsnTMin)) / (pow(((PsnTMax - veg->PsnTMin) / 2.0),2));

	if ((clim->tmin[rstep] < 6) && (DTemp > 0) && (share->GDDTot >= veg->GDDFolEnd))
	{
		//Frost effect
		DTemp = DTemp * (1.0 - ((6.0 - clim->tmin[rstep]) / 6.0) * (share->dayspan / 30.0));  // need to verify zzx
	}

	if (DTemp < 0)
	{
		DTemp = 0;
	}

	share->DVPD = 1.0 - veg->DVPD1 * (pow(share->VPD,veg->DVPD2));
	ten9 = 1000000000.0;

	//Set atmospheric CO2 concentration
	Ca = clim->CO2[rstep];

	//CO2 effect on photosynthesis
	//Leaf internal/external CO2
	CiCaRatio = (-0.075 * veg->FolNCon) + 0.875;
	//Ci at present (350 ppm) CO2
	Ci350 = 350 * CiCaRatio;
	//Ci at RealYear CO2 level

	CiElev = Ca * CiCaRatio;

	//Areal - rate of photosynthesis at a given atmospheric CO2
	//concentration (Ca) relative to that which occurs at CO2 saturation
	Arel350 = 1.22 * ((Ci350 - 68) / (Ci350 + 136));
	ArelElev = 1.22 * ((CiElev - 68) / (CiElev + 136));
	share->DelAmax = 1 + ((ArelElev - Arel350) / Arel350);

	//Calculate CO2 effect on conductance and set slope and intercept for A-gs relationship
	if (site->CO2gsEffect == 1)
	{
		Delgs = share->DelAmax / ((Ca - CiElev) / (350 - Ci350));
		share->DWUE = 1 + (1 - Delgs);					// used for effect on water use efficiency
		gsSlope = (-1.1309 * share->DelAmax) + 1.9762;  //used to determine conductance and then ozone uptake
		gsInt = (0.4656 * share->DelAmax) - 0.9701;
	}
	else
	{
		share->DWUE = 1;
		gsSlope = (-0.6157 * share->DelAmax) + 1.4582;
		gsInt = (0.4974 * share->DelAmax) - 0.9893;
	}




	Amax = (veg->AmaxA + veg->AmaxB * veg->FolNCon) * share->DelAmax;  //nmole CO2/g Fol.sec

	BaseFolResp = veg->BaseFolRespFrac * Amax;
	Amax = Amax * veg->AmaxFrac;
	GrossAmax = Amax + BaseFolResp;
	GrossAmax = (GrossAmax * share->DVPD * DTemp * share->DayLength * 12.0) / ten9;  //g C/g Fol/day, 12 for C atom 

	if (GrossAmax < 0)
	{
		 GrossAmax = 0;
	}

	share->DayResp = ( BaseFolResp * (pow(veg->RespQ10,((share->Tday - veg->PsnTOpt) / 10.0))) * share->DayLength * 12.0) / ten9;
	share->NightResp = (BaseFolResp * (pow(veg->RespQ10,((share->Tnight - veg->PsnTOpt) / 10.0))) * share->NightLength * 12.0) / ten9; //   g C/g Fol/day

	//Initialize ozone effect
	CanopyNetPsnO3 = 0;
	CanopyNetPsnPot = 0;

	//Calculate canopy ozone extinction based on folmass
	O3Prof = 0.6163 + (0.00105 * share->FolMass);


	if (share->FolMass > 0)
	{
		share->CanopyNetPsn = 0;
		share->CanopyGrossPsn = 0;
		share->LAI = 0;
		share->PosCBalMass = share->FolMass;  // posCBalMass: possible C Mass at balance point 
//		O3Effect = 0;
		Layer = 0; 


		for (ix = 1; ix <= 50; ix++)
		{

			i = ix * (share->FolMass / 50.0);
			SLWLayer = veg->SLWmax - (veg->SLWdel * i);//SLW for this layer
			share->LAI = share->LAI + (share->FolMass / 50.0) / SLWLayer;
			Il = clim->par[rstep] * exp(-veg->k * share->LAI );
			LightEff = (1.0 - exp(-Il * log(2.0) / veg->HalfSat));
			LayerGrossPsnRate = GrossAmax * LightEff;					//gC/g Fol/day
			LayerGrossPsn = LayerGrossPsnRate * (share->FolMass / 50.0);   //g C/m2 ground/day
			LayerResp = (share->DayResp + share->NightResp) * (share->FolMass / 50.0);
			LayerNetPsn = LayerGrossPsn - LayerResp;
			if ((LayerNetPsn < 0) && (share->PosCBalMass == share->FolMass))
			{
				share->PosCBalMass = (ix - 1.0) * (share->FolMass / 50.0);

			}


			share->CanopyNetPsn = share->CanopyNetPsn + LayerNetPsn;
			share->CanopyGrossPsn = share->CanopyGrossPsn + LayerGrossPsn;
	    
			//Ozone effects on Net Psn

			if (clim->O3[rstep]>0)
			{
				// Convert netpsn to micromoles for calculating conductance

				netPsnumol = ((LayerNetPsn * 1000000) / (share->DayLength * 12)) / ((share->FolMass / 50) / SLWLayer); //Matlab to C
				//Calculate ozone extinction throughout the canopy
				Layer = Layer + 1;
				RelLayer = Layer / 50.;
				RelO3 = 1 - (RelLayer * O3Prof) * (RelLayer * O3Prof) * (RelLayer * O3Prof); //Matlab to C
				if (RelO3<0) RelO3=0;  // improved by zzx

				//Calculate Conductance (mm/s): Conductance down-regulates with prior O3 effects on Psn
				LayerG = (gsInt + (gsSlope * netPsnumol)) * (1 - share->O3Effect[(int)Layer]);
				//For no downregulation use:    LayerG = gsInt + (gsSlope * netPsnumol);
				if (LayerG < 0)
				{
					LayerG = 0;
				}
				//Calculate cumulative ozone effect for each canopy layer with consideration that previous O3 effects were modified by drought
				share->O3Effect[(int)Layer] = (share->O3Effect[(int)Layer] * share->DroughtO3Frac) + (0.0026 * LayerG * clim->O3[rstep] * RelO3);// 
				if (share->O3Effect[(int)Layer]>1.0)share->O3Effect[(int)Layer]=1.0;
				LayerDO3 = 1 - share->O3Effect[(int)Layer];
			}
			else
			{
				LayerDO3 = 1;
			}

			LayerNetPsnO3 = LayerNetPsn * LayerDO3;
			CanopyNetPsnO3 = CanopyNetPsnO3 + LayerNetPsnO3;




		}// ix=51


		if ((DTemp > 0) && (share->GDDTot > veg->GDDFolEnd) && (clim->doy[rstep] < veg->SenescStart))
		{
			share->PosCBalMassTot = share->PosCBalMassTot + (share->PosCBalMass * share->dayspan);
			share->PosCBalMassIx = share->PosCBalMassIx + share->dayspan;

		}
		if (share->LightEffMin > LightEff)
		{
			share->LightEffMin = LightEff;
		}
	}
	else
	{
		share->PosCBalMass = 0;
		share->CanopyNetPsn = 0;
		share->CanopyGrossPsn = 0;
		share->LAI = 0;
		share->DayResp = 0;
		share->NightResp = 0;
	}

	//Calculate whole-canopy ozone effects before drought
	if ((clim->O3[rstep]>0) && (share->CanopyGrossPsn>0))
	{
		CanopyNetPsnPot = share->CanopyGrossPsn - (share->DayResp * share->FolMass) - (share->NightResp * share->FolMass);
		share->CanopyDO3Pot = CanopyNetPsnO3 / CanopyNetPsnPot;
	}
	else
	{
		share->CanopyDO3Pot = 1;
	}
}




