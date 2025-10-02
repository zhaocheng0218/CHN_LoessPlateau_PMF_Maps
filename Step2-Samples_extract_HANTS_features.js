/*****************************************************************************************************
 * This demo shows the process of using samples to extract HANTS features of Sentinel-2 time series. 
******************************************************************************************************/

var roi = ee.Geometry.Polygon(
        [[[107.46108539811854, 35.50316484800783],
          [107.46108539811854, 34.55181204743065],
          [108.90853412858729, 34.55181204743065],
          [108.90853412858729, 35.50316484800783]]], null, false);
Map.addLayer(roi,{color:"gray"},"roi");
Map.centerObject(roi,11);


////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
////+++ Load Sentinel-2 images.
var start_date = '2020-03-01';
var end_date = '2020-11-01';
var oldBands = ["B2","B3","B4","B8","B11","B12"];
var newBands = ["blue","green","red","nir","swir1","swir2"];

function maskS2clouds(image) {
  var qa = image.select('QA60');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image.updateMask(mask);
}


////+++ Function used for calculate vegetation indices (VIs).
function VIs(img){
  var gcvi = img.expression("(nir/green)-1",{"nir":img.select("nir"),
                                           "green":img.select("green")}).rename("gcvi");
  var ndvi = img.normalizedDifference(["nir","red"]).rename("ndvi");
  var lswi = img.normalizedDifference(["nir","swir1"]).rename("lswi");
  var nmdi = img.expression("(nir-swir1+swir2)/(nir+swir1-swir2)",
                             {"nir":img.select("nir"),
                              "swir1":img.select("swir1"),
                              "swir2":img.select("swir2")
                             }).rename("nmdi");
  var bsi = img.expression("(swir2+red-nir+blue)/(swir2+red+nir-blue)",
                            {"swir2":img.select("swir2"),
                              "red":img.select("red"),
                              "nir":img.select("nir"),
                              "blue":img.select("blue")
                            }).rename("bsi");
  var dbsi = img.normalizedDifference(["swir1","green"]).subtract(ndvi).rename("dbsi");
  var pmli = img.normalizedDifference(["swir1","red"]).rename("pmli");
  
  return img.addBands([gcvi, ndvi, lswi, nmdi, bsi, dbsi, pmli]);
}


var s2Col = ee.ImageCollection("COPERNICUS/S2_SR_HARMONIZED")
              .filterDate(start_date,end_date)
              .filterBounds(roi);
var s2Col = s2Col.map(maskS2clouds)
                .select(oldBands,newBands)
                .map(function(img){return img.divide(10000).copyProperties(img,['system:time_start'])})
                .map(VIs);

////++++Calculate MBPMFI and BPMFI.
var median = s2Col.median().select("blue");  ///++roughly fill gaps

var img_PMS = s2Col.filterDate("2020-03-01","2020-04-01").min().select("blue").unmask(median);
var img_MS =  s2Col.filterDate("2020-05-01","2020-06-01").max().select("blue").unmask(median);
var img_FS =  s2Col.filterDate("2020-07-01","2020-08-01").min().select("blue").unmask(median);

var mbpmfi = img_MS.rename("mbpmfi");
var bpmfi = (img_MS.subtract(img_PMS)).multiply(img_MS.subtract(img_FS))
            .multiply(ee.Number(100)).rename("bpmfi");
////+++++++++++++++++++++++++++++



////+++ Add time band and constant band.
function addVaribles(img){
  var date_img = ee.Date(img.get("system:time_start"));
  var years_img = date_img.difference(ee.Date("2020-01-01"), "year");
  var t_band = ee.Image(years_img).toFloat().rename("t");
  var constant_band = ee.Image(1).toFloat().rename("constant");
  return img.addBands([t_band, constant_band]);
}

var s2Col = s2Col.map(addVaribles);
// print("s2Col:",s2Col.limit(5));
////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
////+++ Harmonic Analysis of Time Series.
////+++ 2 terms, w=1.5 
function addSinCos(img){
  var timeRadians = img.select("t").multiply(2*Math.PI*1.5);
  var timeRadians2 = img.select("t").multiply(4*Math.PI*1.5);
  var cos = timeRadians.cos().rename("cos");
  var sin = timeRadians.sin().rename("sin");
  var cos2 = timeRadians2.cos().rename("cos2");
  var sin2 = timeRadians2.sin().rename("sin2");
  return img.addBands([cos,sin,cos2,sin2]);
}


var s2Col = s2Col.map(addSinCos);
print("s2Col:",s2Col.limit(5));
var harmonicIndependents = ee.List(['constant','sin','cos','sin2','cos2']);
var harmonicDependents = ee.List(['blue','green','red','swir1','swir2',
                                  'gcvi','ndvi','lswi','nmdi','bsi','dbsi','pmli']);

var parameters = harmonicIndependents.add(harmonicDependents).flatten();

var harmonicTrend = s2Col.select(parameters)
                                .reduce(ee.Reducer.linearRegression(harmonicIndependents.size(),
                                                                    harmonicDependents.size()));


var coeffi = harmonicTrend.select("coefficients");
// Map.addLayer(coeffi,{},"coeffi",false);

////+++ coefficient and fit in each band.
var blue_coeffi = coeffi.arraySlice(1,0,1).arrayProject([0]).arrayFlatten([harmonicIndependents])
                        .select(['constant','sin','cos','sin2','cos2'],
                                ['blue_constant','blue_sin','blue_cos','blue_sin2','blue_cos2']);
// Map.addLayer(blue_coeffi,{},"blue_coeffi",false);
var blue_peak = s2Col.select(['constant','sin','cos','sin2','cos2','t'])
                    .map(function(img){
                      return img.addBands(img.select(harmonicIndependents).multiply(blue_coeffi)
                                                            .reduce("sum").rename("blue_fitted"));
                    }).select(["blue_fitted","t"]);
Map.addLayer(blue_peak.select("blue_fitted"),{},"blue_peak",false);
var blue_peak = blue_peak.reduce(ee.Reducer.max(2))
                      .select(["max","max1"],["peak_blue","timing_blue"]);
var green_coeffi = coeffi.arraySlice(1,1,2).arrayProject([0]).arrayFlatten([harmonicIndependents])
                        .select(['constant','sin','cos','sin2','cos2'],
                                ['green_constant','green_sin','green_cos','green_sin2','green_cos2']);
var green_peak = s2Col.select(['constant','sin','cos','sin2','cos2','t'])
                    .map(function(img){
                      return img.addBands(img.select(harmonicIndependents).multiply(green_coeffi)
                                                            .reduce("sum").rename("green_fitted"));
                    }).select(["green_fitted","t"]).reduce(ee.Reducer.max(2))
                      .select(["max","max1"],["peak_green","timing_green"]);
var red_coeffi = coeffi.arraySlice(1,2,3).arrayProject([0]).arrayFlatten([harmonicIndependents])
                        .select(['constant','sin','cos','sin2','cos2'],
                                ['red_constant','red_sin','red_cos','red_sin2','red_cos2']);
var red_peak = s2Col.select(['constant','sin','cos','sin2','cos2','t'])
                    .map(function(img){
                      return img.addBands(img.select(harmonicIndependents).multiply(red_coeffi)
                                                            .reduce("sum").rename("red_fitted"));
                    }).select(["red_fitted","t"]).reduce(ee.Reducer.max(2))
                      .select(["max","max1"],["peak_red","timing_red"]);
var swir1_coeffi = coeffi.arraySlice(1,3,4).arrayProject([0]).arrayFlatten([harmonicIndependents])
                        .select(['constant','sin','cos','sin2','cos2'],
                                ['swir1_constant','swir1_sin','swir1_cos','swir1_sin2','swir1_cos2']);
var swir1_peak = s2Col.select(['constant','sin','cos','sin2','cos2','t'])
                    .map(function(img){
                      return img.addBands(img.select(harmonicIndependents).multiply(swir1_coeffi)
                                                            .reduce("sum").rename("swir1_fitted"));
                    }).select(["swir1_fitted","t"]).reduce(ee.Reducer.max(2))
                      .select(["max","max1"],["peak_swir1","timing_swir1"]);
var swir2_coeffi = coeffi.arraySlice(1,4,5).arrayProject([0]).arrayFlatten([harmonicIndependents])
                        .select(['constant','sin','cos','sin2','cos2'],
                                ['swir2_constant','swir2_sin','swir2_cos','swir2_sin2','swir2_cos2']);
var swir2_peak = s2Col.select(['constant','sin','cos','sin2','cos2','t'])
                    .map(function(img){
                      return img.addBands(img.select(harmonicIndependents).multiply(swir2_coeffi)
                                                            .reduce("sum").rename("swir2_fitted"));
                    }).select(["swir2_fitted","t"]).reduce(ee.Reducer.max(2))
                      .select(["max","max1"],["peak_swir2","timing_swir2"]);

var gcvi_coeffi = coeffi.arraySlice(1,5,6).arrayProject([0]).arrayFlatten([harmonicIndependents])
                        .select(['constant','sin','cos','sin2','cos2'],
                                ['gcvi_constant','gcvi_sin','gcvi_cos','gcvi_sin2','gcvi_cos2']);
var gcvi_peak = s2Col.select(['constant','sin','cos','sin2','cos2','t'])
                    .map(function(img){
                      return img.addBands(img.select(harmonicIndependents).multiply(gcvi_coeffi)
                                                            .reduce("sum").rename("gcvi_fitted"));
                    }).select(["gcvi_fitted","t"]).reduce(ee.Reducer.max(2))
                      .select(["max","max1"],["peak_gcvi","timing_gcvi"]);
var ndvi_coeffi = coeffi.arraySlice(1,6,7).arrayProject([0]).arrayFlatten([harmonicIndependents])
                        .select(['constant','sin','cos','sin2','cos2'],
                                ['ndvi_constant','ndvi_sin','ndvi_cos','ndvi_sin2','ndvi_cos2']);
var ndvi_peak = s2Col.select(['constant','sin','cos','sin2','cos2','t'])
                    .map(function(img){
                      return img.addBands(img.select(harmonicIndependents).multiply(ndvi_coeffi)
                                                            .reduce("sum").rename("ndvi_fitted"));
                    }).select(["ndvi_fitted","t"]).reduce(ee.Reducer.max(2))
                      .select(["max","max1"],["peak_ndvi","timing_ndvi"]);


var lswi_coeffi = coeffi.arraySlice(1,7,8).arrayProject([0]).arrayFlatten([harmonicIndependents])
                        .select(['constant','sin','cos','sin2','cos2'],
                                ['lswi_constant','lswi_sin','lswi_cos','lswi_sin2','lswi_cos2']);
var lswi_peak = s2Col.select(['constant','sin','cos','sin2','cos2','t'])
                    .map(function(img){
                      return img.addBands(img.select(harmonicIndependents).multiply(lswi_coeffi)
                                                            .reduce("sum").rename("lswi_fitted"));
                    }).select(["lswi_fitted","t"]).reduce(ee.Reducer.max(2))
                      .select(["max","max1"],["peak_lswi","timing_lswi"]);
var nmdi_coeffi = coeffi.arraySlice(1,8,9).arrayProject([0]).arrayFlatten([harmonicIndependents])
                        .select(['constant','sin','cos','sin2','cos2'],
                                ['nmdi_constant','nmdi_sin','nmdi_cos','nmdi_sin2','nmdi_cos2']);
var nmdi_peak = s2Col.select(['constant','sin','cos','sin2','cos2','t'])
                    .map(function(img){
                      return img.addBands(img.select(harmonicIndependents).multiply(nmdi_coeffi)
                                                            .reduce("sum").rename("nmdi_fitted"));
                    }).select(["nmdi_fitted","t"]).reduce(ee.Reducer.max(2))
                      .select(["max","max1"],["peak_nmdi","timing_nmdi"]);

var bsi_coeffi = coeffi.arraySlice(1,9,10).arrayProject([0]).arrayFlatten([harmonicIndependents])
                        .select(['constant','sin','cos','sin2','cos2'],
                                ['bsi_constant','bsi_sin','bsi_cos','bsi_sin2','bsi_cos2']);
var bsi_peak = s2Col.select(['constant','sin','cos','sin2','cos2','t'])
                    .map(function(img){
                      return img.addBands(img.select(harmonicIndependents).multiply(bsi_coeffi)
                                                            .reduce("sum").rename("bsi_fitted"));
                    }).select(["bsi_fitted","t"]).reduce(ee.Reducer.max(2))
                      .select(["max","max1"],["peak_bsi","timing_bsi"]);
var dbsi_coeffi = coeffi.arraySlice(1,10,11).arrayProject([0]).arrayFlatten([harmonicIndependents])
                        .select(['constant','sin','cos','sin2','cos2'],
                                ['dbsi_constant','dbsi_sin','dbsi_cos','dbsi_sin2','dbsi_cos2']);
var dbsi_peak = s2Col.select(['constant','sin','cos','sin2','cos2','t'])
                    .map(function(img){
                      return img.addBands(img.select(harmonicIndependents).multiply(dbsi_coeffi)
                                                            .reduce("sum").rename("dbsi_fitted"));
                    }).select(["dbsi_fitted","t"]).reduce(ee.Reducer.max(2))
                      .select(["max","max1"],["peak_dbsi","timing_dbsi"]);

var pmli_coeffi = coeffi.arraySlice(1,11,12).arrayProject([0]).arrayFlatten([harmonicIndependents])
                        .select(['constant','sin','cos','sin2','cos2'],
                                ['pmli_constant','pmli_sin','pmli_cos','pmli_sin2','pmli_cos2']);
var pmli_peak = s2Col.select(['constant','sin','cos','sin2','cos2','t'])
                    .map(function(img){
                      return img.addBands(img.select(harmonicIndependents).multiply(pmli_coeffi)
                                                            .reduce("sum").rename("pmli_fitted"));
                    }).select(["pmli_fitted","t"]).reduce(ee.Reducer.max(2))
                      .select(["max","max1"],["peak_pmli","timing_pmli"]);
                       
var data = blue_coeffi.addBands(green_coeffi).addBands(red_coeffi).addBands(swir1_coeffi).addBands(swir2_coeffi)
                      .addBands(gcvi_coeffi).addBands(ndvi_coeffi).addBands(lswi_coeffi).addBands(nmdi_coeffi)
                      .addBands(bsi_coeffi).addBands(dbsi_coeffi).addBands(pmli_coeffi)
                      .addBands(green_peak).addBands(red_peak).addBands(swir1_peak).addBands(swir2_peak)
                      .addBands(gcvi_peak).addBands(ndvi_peak).addBands(lswi_peak).addBands(nmdi_peak)
                      .addBands(bsi_peak).addBands(dbsi_peak).addBands(pmli_peak).addBands(blue_peak);
var data = data.addBands([mbpmfi,bpmfi]);
print("data:",data);
Map.addLayer(data,{},"data",false);
////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


////+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
////+++ Use samples to extract features.
var raw_pmfSamples = ee.FeatureCollection("users/my-work/LoessPlateau_PFM/GitHub_openAccess/pmf_Samples")
                    .map(function(fea){return ee.Feature(fea.geometry())})
                    .map(function(fea){return fea.set("landcover",1)});
var raw_nopmfSamples = ee.FeatureCollection("users/my-work/LoessPlateau_PFM/GitHub_openAccess/nopmf_Samples")
                    .map(function(fea){return ee.Feature(fea.geometry())})
                    .map(function(fea){return fea.set("landcover",0)});
                    
Map.addLayer(raw_pmfSamples,{color:"red"},"raw_pmfSamples");
Map.addLayer(raw_nopmfSamples,{color:"blue"},"raw_nopmfSamples");


var pmfSamples = data.sampleRegions({
                      collection:raw_pmfSamples,
                      properties:["landcover"],
                      scale:10,
                      tileScale:16,
                      geometries:true
                    });
var nopmfSamples = data.sampleRegions({
                      collection:raw_nopmfSamples,
                      properties:["landcover"],
                      scale:10,
                      tileScale:16,
                      geometries:true
                    });

print("pmfSamples",pmfSamples.limit(10));
print("nopmfSamples",nopmfSamples.limit(10));

Export.table.toAsset({
  collection: pmfSamples,
  description:"pmf_Samples_HANTS"
});

Export.table.toAsset({
  collection: nopmfSamples,
  description:"nopmf_Samples_HANTS"
});


