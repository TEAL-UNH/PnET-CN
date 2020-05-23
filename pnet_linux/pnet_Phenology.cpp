#include "pnet_model.h"

void pnet_model::Phenology(veg_struct* veg, clim_struct* clim, int rstep, share_struct* share, int GrowthPhase)
{
	//
	// Phenlology calculations for the PnET ecosystem model.
	//

	double OldFolMass, GDDFolEff, delGDDFolEff, FolMassNew,GDD;


	if (GrowthPhase == 1)
	{

		share->dayspan = (double)getdays(clim->year[rstep],clim->doy[rstep]);
//		share->dayspan = 30;   // matlab version test

		GDD = share->Tave * share->dayspan;

		share->TaveYr += GDD/365.0 ;
		share->PARYr += clim->par[rstep]* share->dayspan/365.0 ;

		if (GDD < 0 )  GDD = 0;
		share->GDDTot = share->GDDTot + GDD;

		if ((share->GDDTot >= veg->GDDFolStart) && (clim->doy[rstep] < veg->SenescStart))
		{
			OldFolMass = share->FolMass;  
			GDDFolEff = (share->GDDTot - veg->GDDFolStart) / (veg->GDDFolEnd - veg->GDDFolStart);

			if (GDDFolEff < 0)GDDFolEff = 0;
			if (GDDFolEff > 1)GDDFolEff = 1;

			delGDDFolEff = GDDFolEff - share->OldGDDFolEff;
			share->FolMass = share->FolMass + (share->BudC * delGDDFolEff) / veg->CFracBiomass;
			share->FolProdCMo = (share->FolMass - OldFolMass) * veg->CFracBiomass;
			share->FolGRespMo = share->FolProdCMo * veg->GRespFrac;
			share->OldGDDFolEff = GDDFolEff;
		}
		else
		{
			share->FolProdCMo = 0; 
			share->FolGRespMo = 0;
		}
	}
	else
	{
		share->FolLitM = 0;
		if ((share->PosCBalMass < share->FolMass) && (clim->doy[rstep] > veg->SenescStart))
		{


			if ((share->PosCBalMass) > (veg->FolMassMin))
			{
				FolMassNew = share->PosCBalMass;
			}
			else
			{
				FolMassNew = veg->FolMassMin;
			}


			if (FolMassNew == 0)
			{
				share->LAI = 0;
			}
			else if (FolMassNew < share->FolMass)
			{
				share->LAI = share->LAI * (FolMassNew / share->FolMass);
			}
			if (FolMassNew < share->FolMass)
			{
				share->FolLitM = share->FolMass - FolMassNew;
			}
			share->FolMass = FolMassNew;
		}
	}


}
