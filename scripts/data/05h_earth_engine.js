// =============================================================================
// 05H_earth_engine_covariates.js
// =============================================================================
// Run in Google Earth Engine Code Editor after uploading the canonical GeoJSON.
// Exports source-level district summaries; anomalies/lags are constructed in R
// by 05L using fixed pre-study baselines where possible.
//
// Exports to Google Drive folder: cabo_delgado_covariates
//   05H_viirs_monthly_2012_2024.csv
//   05H_modis_vi_monthly_2001_2024.csv
//   05H_era5land_monthly_2017_2024.csv
//   05H_dynamicworld_annual_2016_2024.csv
//   05H_hansen_annual_2016_2024.csv
// =============================================================================

// EDIT THIS after uploading 05G GeoJSON as an Earth Engine table asset.
var DISTRICT_ASSET = 'projects/YOUR_PROJECT/assets/cabo_delgado_admin2_canonical';
var DRIVE_FOLDER = 'cabo_delgado_covariates';

var districts = ee.FeatureCollection(DISTRICT_ASSET)
  .select(['district_id', 'canonical_name', 'district_index']);
print('District count (must be 17):', districts.size());

function monthStarts(start, endExclusive) {
  start = ee.Date(start);
  endExclusive = ee.Date(endExclusive);
  var n = endExclusive.difference(start, 'month').toInt();
  return ee.List.sequence(0, n.subtract(1)).map(function(k) {
    return start.advance(ee.Number(k), 'month');
  });
}

function yearList(startYear, endYear) {
  return ee.List.sequence(startYear, endYear);
}

function reduceMeanByDistrict(image, date, scale) {
  date = ee.Date(date);
  return districts.map(function(f) {
    var stats = image.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: f.geometry(),
      scale: scale,
      maxPixels: 1e10,
      bestEffort: false,
      tileScale: 4
    });
    return ee.Feature(null, stats)
      .set('district_id', f.get('district_id'))
      .set('canonical_name', f.get('canonical_name'))
      .set('district_index', f.get('district_index'))
      .set('year_month', date.format('YYYY-MM-01'))
      .set('year', date.get('year'))
      .set('month', date.get('month'));
  });
}

function flattenMappedFeatureCollections(listOfFCs) {
  return ee.FeatureCollection(listOfFCs).flatten();
}

// -----------------------------------------------------------------------------
// A. VIIRS monthly nighttime lights, 2012-04 through 2024-12
// Keep cf_cvg because avg_rad == 0 can mean inadequate usable observations.
// -----------------------------------------------------------------------------
var viirs = ee.ImageCollection('NOAA/VIIRS/DNB/MONTHLY_V1/VCMCFG');
var viirsMonths = monthStarts('2012-04-01', '2025-01-01');
var viirsFC = flattenMappedFeatureCollections(viirsMonths.map(function(ds) {
  ds = ee.Date(ds);
  var de = ds.advance(1, 'month');
  var one = viirs.filterDate(ds, de).mean();
  var cov = one.select('cf_cvg').rename('ntl_cf_cvg');
  var rad = one.select('avg_rad')
    .updateMask(cov.gt(0))
    .rename('ntl_avg_rad');
  var image = rad.addBands(cov);
  return reduceMeanByDistrict(image, ds, 500);
}));

Export.table.toDrive({
  collection: viirsFC,
  description: '05H_viirs_monthly_2012_2024',
  folder: DRIVE_FOLDER,
  fileNamePrefix: '05H_viirs_monthly_2012_2024',
  fileFormat: 'CSV'
});

// -----------------------------------------------------------------------------
// B. MODIS monthly NDVI/EVI, 2001-01 through 2024-12.
// Fixed pre-study 2001-2016 climatology will be computed by 05L.
// QA mask: SummaryQA 0 (good) or 1 (marginal); scale 0.0001.
// -----------------------------------------------------------------------------
var modis = ee.ImageCollection('MODIS/061/MOD13A3');
var modisMonths = monthStarts('2001-01-01', '2025-01-01');
var modisFC = flattenMappedFeatureCollections(modisMonths.map(function(ds) {
  ds = ee.Date(ds);
  var de = ds.advance(1, 'month');
  var c = modis.filterDate(ds, de);
  var img = ee.Image(c.first());
  var qa = img.select('SummaryQA');
  var good = qa.lte(1);
  var vi = img.select(['NDVI', 'EVI'])
    .multiply(0.0001)
    .updateMask(good)
    .rename(['ndvi', 'evi']);
  return reduceMeanByDistrict(vi, ds, 1000);
}));

Export.table.toDrive({
  collection: modisFC,
  description: '05H_modis_vi_monthly_2001_2024',
  folder: DRIVE_FOLDER,
  fileNamePrefix: '05H_modis_vi_monthly_2001_2024',
  fileFormat: 'CSV'
});

// -----------------------------------------------------------------------------
// C. ERA5-Land monthly, secondary climate sensitivity.
// temperature_2m K -> C; volumetric soil water layer 1 retained in m3/m3.
// -----------------------------------------------------------------------------
var era5 = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR');
var eraMonths = monthStarts('2017-01-01', '2025-01-01');
var eraFC = flattenMappedFeatureCollections(eraMonths.map(function(ds) {
  ds = ee.Date(ds);
  var de = ds.advance(1, 'month');
  var img = era5.filterDate(ds, de).first();
  img = ee.Image(img);
  var temp = img.select('temperature_2m').subtract(273.15).rename('era5_temperature_c');
  var soil = img.select('volumetric_soil_water_layer_1').rename('era5_soil_water_l1');
  return reduceMeanByDistrict(temp.addBands(soil), ds, 11132);
}));

Export.table.toDrive({
  collection: eraFC,
  description: '05H_era5land_monthly_2017_2024',
  folder: DRIVE_FOLDER,
  fileNamePrefix: '05H_era5land_monthly_2017_2024',
  fileFormat: 'CSV'
});

// -----------------------------------------------------------------------------
// D. Dynamic World annual land-cover shares, 2016-2024 (2016 supports a 1y lag for 2017).
// Annual modal label. Label codes: water=0, trees=1, grass=2,
// flooded_vegetation=3, crops=4, shrub_and_scrub=5, built=6, bare=7, snow=8.
// Zonal mean of 0/1 band = area share approximately (equal-area caveat is tiny
// at district scale; pixel area weighting can be added as a sensitivity).
// -----------------------------------------------------------------------------
var dw = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1');
var dwFC = flattenMappedFeatureCollections(yearList(2016, 2024).map(function(y) {
  y = ee.Number(y);
  var ds = ee.Date.fromYMD(y, 1, 1);
  var de = ds.advance(1, 'year');
  var label = dw.filterDate(ds, de).select('label').mode();
  var image = label.eq(1).rename('forest_share')
    .addBands(label.eq(4).rename('cropland_share'))
    .addBands(label.eq(6).rename('built_share'))
    .addBands(label.eq(0).rename('water_share'))
    .addBands(label.eq(5).rename('shrub_share'))
    .addBands(label.eq(7).rename('bare_share'));
  return reduceMeanByDistrict(image, ds, 10);
}));

Export.table.toDrive({
  collection: dwFC,
  description: '05H_dynamicworld_annual_2016_2024',
  folder: DRIVE_FOLDER,
  fileNamePrefix: '05H_dynamicworld_annual_2016_2024',
  fileFormat: 'CSV'
});

// -----------------------------------------------------------------------------
// E. Hansen Global Forest Change v1.13, 2016-2024 (2016 supports a 1y lag).
// treecover2000 = baseline percent canopy cover. lossyear 1..25 corresponds to
// 2001..2025. Annual output includes annual and cumulative forest-loss share.
// -----------------------------------------------------------------------------
var hansen = ee.Image('UMD/hansen/global_forest_change_2025_v1_13');
var tree = hansen.select('treecover2000').rename('treecover2000_pct');
var lossYear = hansen.select('lossyear');
var hansenFC = flattenMappedFeatureCollections(yearList(2016, 2024).map(function(y) {
  y = ee.Number(y);
  var code = y.subtract(2000);
  var ds = ee.Date.fromYMD(y, 1, 1);
  var annualLoss = lossYear.eq(code).rename('forest_loss_share_annual');
  var cumulativeLoss = lossYear.gt(0).and(lossYear.lte(code)).rename('forest_loss_share_cumulative');
  var image = tree.addBands(annualLoss).addBands(cumulativeLoss);
  return reduceMeanByDistrict(image, ds, 30);
}));

Export.table.toDrive({
  collection: hansenFC,
  description: '05H_hansen_annual_2016_2024',
  folder: DRIVE_FOLDER,
  fileNamePrefix: '05H_hansen_annual_2016_2024',
  fileFormat: 'CSV'
});

print('Create/Run all five export tasks from the Tasks tab.');
