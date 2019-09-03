IOdictionary sparkProperties
(
	IOobject
    (
		"sparkProperties",
		runTime.constant(),
		runTime,
		IOobject::MUST_READ_IF_MODIFIED,
		IOobject::NO_WRITE
	)
);

volScalarField sparkPVs
(
	IOobject
    (
     	"sparkingPVs",
         mesh.time().timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
    ),
    0*meanPVs
);

Switch sparkHotSpot(sparkProperties.lookup("sparkHotSpot"));

dimensionedScalar sparkBeg("sparkBeg", dimTime, 0.0);
dimensionedScalar sparkDuration("sparkDuration", dimTime, 0.0);
dimensionedScalar sparkingPVs("sparkingPVs", dimensionSet(1, -3, -1, 0, 0, 0, 0), 0.0);
dimensionedVector sparkCenter = 0*mesh.C()[0];
dimensionedScalar sparkRadi("sparkRadi",dimLength, 0.0);
dimensionedScalar sparkL("sparkL",dimLength, 0.0);
dimensionedScalar sparkW("sparkW",dimLength, 0.0);
dimensionedScalar sparkH("sparkH",dimLength, 0.0);
dimensionedScalar minSparkingF("minSparkingF",dimless, 0.0);
dimensionedScalar maxSparkingF("maxSparkingF",dimless, 0.0);
dimensionedScalar sparkingT("maxSparkingT",dimTemperature, 0.0);

bool spherSpark = false;
bool boxSpark   = false;
 
if(sparkHotSpot)
{
sparkProperties.lookup("sparkBeg"  )         >> sparkBeg;
sparkProperties.lookup("sparkDuration"  )    >> sparkDuration;
sparkProperties.lookup("sparkingPVs")   >> sparkingPVs;
sparkProperties.lookup("sparkCenter")        >> sparkCenter;
sparkProperties.lookup("spherSpark")  >> spherSpark;
sparkProperties.lookup("boxSpark")  >> boxSpark;
sparkProperties.lookup("maxSparkingF")  >> maxSparkingF;
sparkProperties.lookup("minSparkingF")  >> minSparkingF;
sparkProperties.lookup("sparkingT")  >> sparkingT;
if(spherSpark)
{
	sparkProperties.lookup("sparkRadi")  >> sparkRadi;
}
if(boxSpark)
{
sparkProperties.lookup("sparkL")  >> sparkL;
sparkProperties.lookup("sparkW")  >> sparkW;
sparkProperties.lookup("sparkH")  >> sparkH;
}

}
