{
     wordList species(combustion->lookup("species"));
     if(runTime.outputTime())
     {
	if(species.size()>0)
        {
	  forAll(species,i)
          {
	     label k = thermo.composition().species()[species[i]];
	     thermo.composition().Y()[k].write();
	  }
	}
	else
	{
	  PtrList<volScalarField>& Y(thermo.composition().Y());
          forAll(Y,i)
          {
            thermo.composition().Y()[i].write();
          }
	}
     }
}
