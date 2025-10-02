/*****************************************************************************************************
 * There is a demo illustrating the automatic training sample generation method.
******************************************************************************************************/

////++Select a typical region for mapping
var roi = ee.Geometry.Polygon(
        [[[107.46108539811854, 35.50316484800783],
          [107.46108539811854, 34.55181204743065],
          [108.90853412858729, 34.55181204743065],
          [108.90853412858729, 35.50316484800783]]], null, false);
Map.addLayer(roi,{color:"red"},"roi",false);
Map.centerObject(roi,12);

////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
////++ Overlap two land use/cover products to exclude no-PMF pixels.
////++ Dataset 1: GLC-FCS30D (Zhang et al., 2024)
////++ Dataset 2: CACD (Tu et al., 2024)
////++ Dataset 3: CLCD (Huang et al., 2021)
var cropland = ee.Image("users/my-work/LoessPlateau_PFM/GitHub_openAccess/cropland")
Map.addLayer(cropland.selfMask().clip(roi),{palette:"#d774ff"},"cropland",false);
////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


/////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/////++ Load data.
var start_date = '2020-02-01';
var end_date = '2020-10-01';
var oeel=require('users/OEEL/lib:loadAll');
var oldBands = ["B2","B3","B4","B5","B6","B7","B8","B8A","B11","B12"];
var newBands = ["blue","green","red","re1","re2","re3","nir","re4","swir1","swir2"];


////+++Sentinel-2 data.
var s2Col = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
              .filterDate(start_date,end_date)
              .filterBounds(roi)
              .map(function(img){return img.clip(roi)});
    s2Col = oeel.Algorithms.Sentinel2.cloudfree(50,s2Col)
                .select(oldBands,newBands)
                .map(function(img){return img.divide(10000)
                                  .copyProperties(img,["system:time_start","system:index"]); })
                .select(["red","green","blue"]);

////++ High-quality Sentinel-2 images during pre-mulching stage (PMS)
var timing_PMS = ee.Date("2020-03-15").millis();
var img_PMS = s2Col.filterDate("2020-03-15","2020-03-20").min()
                  .set("system:time_start",timing_PMS);
Map.addLayer(img_PMS,{bands:["red","green","blue"],min:0.05,max:0.4},"img_PMS",false);

////++ High-quality Sentinel-2 images during mulching stage (MS)
var timing_max = ee.Date("2020-05-01").millis();
var img_max = s2Col.filter( ee.Filter.or(ee.Filter.date("2020-05-15","2020-05-20"),
                                         ee.Filter.date("2020-05-25","2020-05-30")) )
                  .max()
                  .set("system:time_start",timing_max);
Map.addLayer(img_max,{bands:["red","green","blue"],min:0.017,max:0.27},"img_max");

////++ High-quality Sentinel-2 images during flourishing stage (FS)
var timing_FS = ee.Date("2020-08-15").millis();
var img_FS = s2Col.filterDate("2020-07-01","2020-07-10").min()
                  .set("system:time_start",timing_FS);
Map.addLayer(img_FS,{bands:["red","green","blue"],min:0.05,max:0.4},"img_FS",false);

////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
////++ Calculating plastic-mulched farmland indices.
////+++ MBPMFI = img_max
var MBPMFI = img_max.select("blue").updateMask(cropland).rename("MBPMFI");
////+++ BPMFI = (img_max-img_PMS)*(img_max-img_FS)
var BPMFI = (img_max.subtract(img_PMS)).multiply(img_max.subtract(img_FS))
            .multiply(ee.Number(100)).select("blue")
            .updateMask(cropland).rename("BPMFI");
////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
////++ Threshold-based method for automatically generating traning samples.
var pmf1 = MBPMFI.gte(0.16).eq(1).rename("pmf1");
var pmf2 = BPMFI.gte(0.40).eq(1).rename("pmf2");
var pmf = pmf1.add(pmf2).eq(2);
var pmf = pmf.updateMask(pmf).rename("pmf");
Map.addLayer(pmf,{palette:["red"],min:0,max:1},"pmf");

var no_pmf1 = MBPMFI.lt(0.16).eq(1).rename("no_pmf1");
var no_pmf2 = BPMFI.lt(0.40).eq(1).rename("no_pmf2");
var no_pmf = no_pmf1.add(no_pmf2).eq(2);
var no_pmf = no_pmf.updateMask(no_pmf).rename("nopmf");
Map.addLayer(no_pmf,{palette:["#03fff8"],min:0,max:1},"no_pmf");

////+++ Neighborhood filter.
var kernel = ee.Kernel.square(1,"pixels",false);
var kernelArea = ee.Number(3).pow(2);

//++PMF neighborhood image
var pmf_Neighbor = pmf.unmask(0).clip(roi).convolve(kernel).rename("pmfNei");

//++no-PMF neighborhood image
var nopmf_Neighbor = no_pmf.unmask(0).clip(roi).convolve(kernel).rename("nopmfNei");

//++Randomly stratified sampling
//++PMF
var pmf_Samples = pmf.stratifiedSample({
  numPoints:500,
  region:roi,
  scale:10,
  tileScale:16,
  geometries:true
});
var pmf_Samples = pmf_Neighbor.sampleRegions({
  collection:pmf_Samples,
  properties:["pmf"],
  scale:10,
  tileScale:16,
  geometries:true
});
var pmf_Samples = pmf_Samples.filter(ee.Filter.eq("pmfNei",9));
print("pmf_Samples",pmf_Samples.limit(10));


//++no-PMF
var nopmf_Samples = no_pmf.stratifiedSample({
  numPoints:500,
  region:roi,
  scale:10,
  tileScale:16,
  geometries:true
});
var nopmf_Samples = nopmf_Neighbor.sampleRegions({
  collection:nopmf_Samples,
  properties:["nopmf"],
  scale:10,
  tileScale:16,
  geometries:true
});
var nopmf_Samples = nopmf_Samples.filter(ee.Filter.eq("nopmfNei",9));
print("nopmf_Samples",nopmf_Samples.limit(10));
////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
////++ Export traning samples to Assets.
Export.table.toAsset({
  collection: pmf_Samples,
  description:"pmf_Samples"
});

Export.table.toAsset({
  collection: nopmf_Samples,
  description:"nopmf_Samples"
});

