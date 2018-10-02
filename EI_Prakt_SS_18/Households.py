import arcpy
import arcpy as env
import json
import requests
import numpy as np
import pandas as pd
import numpy.lib.recfunctions as rfn
import math

def extract_buildings(building_data, street_data, building_usage):
    # Datapath for global usage in this function
    # Create layer of buildings
    arcpy.MakeFeatureLayer_management(building_data, 'buildings_data_layer')

    # Create layer of streets
    arcpy.MakeFeatureLayer_management(street_data, 'street_layer')

    # Select all residential buildings
    if building_usage == 'r':
        arcpy.SelectLayerByAttribute_management("buildings_data_layer", "NEW_SELECTION",
                                                "building = 'apartments' "
                                                "or building = 'dormitory' "
                                                "or building = 'house' "
                                                "or building = 'hut' "
                                                "or building = 'residential' ")

    # Select all non residential buildings
    elif building_usage == 'nr':
        arcpy.SelectLayerByAttribute_management("buildings_data_layer", "NEW_SELECTION",
                                                "building = 'Cinema' "
                                                "or building = 'commercial' "
                                                "or building = 'construction' "
                                                "or building = 'hospital' "
                                                "or building = 'hotel' "
                                                "or building = 'museum' "
                                                "or building = 'office' "
                                                
                                                "or building = 'public' "
                                                "OR building = 'restaurant' "
                                                "OR building = 'retail' "
                                                "or building = 'school' "
                                                "or building = 'train_station' "
                                                "or building = 'university' ")
        # or building = yes and amenity = ....

    # Select all buildings
    else:
        arcpy.SelectLayerByAttribute_management("buildings_data_layer", "NEW_SELECTION",
                                                "building = 'apartments' "
                                                "or building = 'cathedral' "
                                                "or building = 'church' "
                                                "or building = 'Cinema' "
                                                "or building = 'collapsed' "
                                                "or building = 'commercial' "
                                                "or building = 'construction' "
                                                "or building = 'dormitory' "
                                                "or building = 'garage' "
                                                "or building = 'garages' "
                                                "or building = 'hospital' "
                                                "or building = 'hotel' "
                                                "or building = 'house' "
                                                "or building = 'hut' "
                                                "or building = 'museum' "
                                                "or building = 'office' "
                                                "or building = 'public' "
                                                "or building = 'residential' "
                                                "OR building = 'restaurant' "
                                                "OR building = 'retail' "
                                                "OR building = 'roof' "
                                                "or building = 'school' "
                                                "or building = 'terrace' "
                                                "or building = 'train_station' "
                                                "or building = 'transportation' "
                                                "or building = 'university' "
                                                "or building = 'yes'")

    # Export buildings and streets as feature classes
    arcpy.FeatureClassToFeatureClass_conversion('buildings_data_layer', r'C:\...\Output', 'buildings_feature_class')
    arcpy.FeatureClassToFeatureClass_conversion('street_layer', r'C:\...\Output', 'street_feature_class')

    # Add geodesic area in m^2
    arcpy.AddGeometryAttributes_management("buildings_feature_class.shp", "AREA_GEODESIC", "METERS", "SQUARE_METERS")

    #building_lvls = r'C:\...\Output\building_lvl'
    #arcpy.JoinField_management("buildings_feature_class.shp", "AREA_GEODESIC", building_lvls, "building:levels")

    # Add  field number of estimated households per building
    arcpy.AddField_management("buildings_feature_class.shp", "households", "SHORT")

    # Save buildings as centroid
    arcpy.FeatureToPoint_management("buildings_feature_class.shp", r"C:\...\Output\buildings_feature_class_pt", "CENTROID")

    # near table map centroids to closest street segment
    arcpy.GenerateNearTable_analysis("buildings_feature_class_pt.shp", "street_feature_class.shp")
    print "Generated requested data"


def bbox_from_extent(feature_class):
    extent = arcpy.Describe(feature_class).extent
    minLon = extent.XMin
    minLat = extent.YMin
    maxLon = extent.XMax
    maxLat = extent.YMax
    return minLat, minLon, maxLat, maxLon


# This function extracts building levels in a given bounding box secstor
def get_building_levels():
    # Datapath for global usage in this function
    global data_test_osm

    # call function to create bbox from extent_fc coordinates
    bbox_west, bbox_south, bbox_east, bbox_north = bbox_from_extent(data_test_osm)

    # bounding box frames are adapted to the extent of the current shapefile
    arcpy.env.extent = "%s %s %s %s" % (bbox_west, bbox_south, bbox_east, bbox_north)
    print("Current bbox set to %s %s %s %s " % (bbox_west, bbox_south, bbox_east, bbox_north))

    # Parameters for GET Request to overpass Turbo for:
    # all buildings in both open and closed ways containing information about their level height
    # within the coordinates of the specified bounding box
    # remove redundancies
    # TODO: recheck OSM request
    osmrequest = {'data': '[out:json][timeout:25];'
                          '('
                          'way["building:levels"]'
                          '(%s, %s, %s, %s);'
                          '(._;>;);'
                          ');'
                          'out;'
                          % (bbox_west, bbox_south, bbox_east, bbox_north)}
    # URL for GET Request
    osmurl = 'http://overpass-api.de/api/interpreter'

    # GET Request to overpass API
    osm = requests.get(osmurl, params=osmrequest)

    # Extract data from JSON dictionary
    osmdata = osm.json()
    #print(osmdata)
    osmdata = osmdata['elements']
    #print(osmdata)
    # rearrange osmdata so only relevant attributes can be extracted as dataframe
    for i in osmdata:
        if 'tags' in i:
            for k, v in i['tags'].iteritems():
                i[k] = v
            del i['tags']

    # Create a dataframe from neatly sorted JSON dictionary containing all the building:level information
    # TODO : 'id'
    osmdataframe = pd.DataFrame(osmdata, columns=['id', 'building:levels'])

    # coerce dtype 'O' object values to np array native float64 type
    osmdataframe['building:levels'] = pd.to_numeric(osmdataframe['building:levels'], errors='coerce')

    # remove  all rows containing NaNs (no building level information)
    osmdataframe = osmdataframe[~np.isnan(osmdataframe).any(axis=1)]

    # Create a sorted np array
    sorted_osmNParray = np.array(osmdataframe.to_records())

    print(sorted_osmNParray.dtype)

    # reconvert to np array so it can be passed to NumpyArrayToTable method
    osmNParray = np.array(sorted_osmNParray, np.dtype([('index', '<i8'), ('id', '<i8'), ('building:levels', '<f8')]))
    print(repr(osmNParray))

    # rename id field so it can be used for a later join
    #osmNParray_ren = rfn.rename_fields(osmNParray, {'id': 'OSMID'})

    # Convert NP Array to table
    arcpy.da.NumPyArrayToTable(osmNParray, r'C:\Users\Andi\Documents\ArcGIS 10.3.1\NumberOfHouseholds\NumberOfHouseholds.gdb\out_table')
    #arcpy.DeleteField_management(r'C:\Users\Andi\Documents\ArcGIS 10.3.1\NumberOfHouseholds\NumberOfHouseholds.gdb\out_table', "index")

    # Join Operation
    # TODO: JOINT OPERATION
    # in data: feature class you want to use
    arcpy.JoinField_management("buildings_feature_class.shp", "OSMID", r'C:\Users\Andi\Documents\ArcGIS 10.3.1\NumberOfHouseholds\NumberOfHouseholds.gdb\out_table', "id", ["building_levels"])
    # Calculation of households
    # avg population per household = population_acrage_density * acrage * building levels = inhabitants per building
    #             row[2] = [0.045 * (row[0] * (row[1]+1))]
    # Rundung in Household Feld, so dass keine strangen Berechnungen raus kommen
    # 10 signs per field maximum
    with arcpy.da.UpdateCursor("buildings_feature_class.shp", ["AREA_GEO", "building_l", "households"]) as cursor:
        for row in cursor:
            # [pop_dnsty_const * AREA_GEO * (building_levels + 1)] / #avg inhabitants per household
            # maybe even a get request to determine the avg population density
            row[2] = (0.045 * row[0] * (row[1]+1))/1.78
            cursor.updateRow(row)
    print "EXIT"


if __name__ == '__main__':
    # set workspace
    outpath = r'C:\...'
    arcpy.env.overwriteOutput = True
    
    # Output
    arcpy.env.workspace = outpath

    # Data containing building information
    data_test_osm = r'C:\...'
    roads_cleared = r'C:\...'

    # get user input
    building_type = raw_input("Please enter \'r\' if you would like to search for residential buildings. \n "
                              "Enter \'nr\' if you would like to search for non residential (i.e. commercial and civic)"
                              "buildings. \n"
                              "Enter \'all\' if you would like to search for all buildings. \n")

    # Check user input for validity
    if building_type != 'r' and building_type != 'nr' and building_type != 'all':
        print "invalid input"
    else:
        # Call function
        print("function call in progress...")
        extract_buildings(data_test_osm, roads_cleared, building_type)
        get_building_levels()
