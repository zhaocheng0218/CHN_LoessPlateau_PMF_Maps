var roi = ee.FeatureCollection("users/my-work/LoessPlateau_PFM/LoessPlateau")
            .first().geometry();
var city = ee.FeatureCollection("users/my-work/LoessPlateau_PFM/LP_city_nominal");

var image = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Baiyin2020"),
    image2 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Baotou2020"),
    image3 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Bayannur2020"),
    image4 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Changzhi2020"),
    image5 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Datong2020"),
    image6 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Dingxi2020"),
    image7 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Guanzhong2020"),
    image8 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Guyuan2020"),
    image9 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Henan2020"),
    image10 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Hohhot2020"),
    image11 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Jincheng2020"),
    image12 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Jinzhong2020"),
    image13 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Lanzhou2020"),
    image14 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Linfen2020"),
    image15 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Linxia2020"),
    image16 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Lvliang2020"),
    image17 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_NID2020"),
    image18 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Odors2020"),
    image19 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Pingliang2020"),
    image20 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Qinghai2020"),
    image21 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Qingyang2020"),
    image22 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Shuozhou2020"),
    image23 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Taiyuan2020"),
    image24 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Tianshui2020"),
    image25 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Ulanqab2020"),
    image26 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Wuzhong2020"),
    image27 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Xinzhou2020"),
    image28 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_YanAn2020"),
    image29 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Yangquan2020"),
    image30 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Yulin2020"),
    image31 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Yuncheng2020"),
    image32 = ee.Image("users/my-work/LoessPlateau_PFM/PFM_results_2020/pfm_Zhongwei2020");

var img = ee.ImageCollection([image,image2,image3,image4,image5,image6,image7,image8,image9,image10,
                              image11,image12,image13,image14,image15,image16,image17,image18,image19,image20,
                              image21,image22,image23,image24,image25,image26,image27,image28,image29,image30,
                              image31,image32]).mosaic();
Map.addLayer(img,{palette:"yellow"},"img",false);

var cityList = city.toList(32);
print(ee.Feature(cityList.get(1)).get("City"))

for (var i=0;i<32;i++){
  var temp_roi = cityList.get(i);
  var temp_name = ee.String(ee.Feature(temp_roi).get("City")).getInfo();
  var temp_geo = ee.Feature(temp_roi).geometry();
  var temp_img = img.clip(temp_geo);
  
  Export.image.toDrive({
  image:temp_img,
  description:temp_name+"_pfm_Results2020",
  region:temp_geo,
  scale:10,
  crs:"EPSG:4326",
  maxPixels:1e13
});

}



