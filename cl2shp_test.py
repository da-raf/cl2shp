import cl2shp
import unittest

import os.path

test_data_file = os.path.join('test_data', 'test_data_file.osm')
rectangle = cl2shp.Rect( [-55.0,-35.0], [-53.5,-33.8] )
num_threads = 4


def getWaysContainingNodeId(ways, node_id):
    result_ways = []
    for way in ways:
        if node_id in way.nodes():
            result_ways.append(way)
    return result_ways


class cl2shpTester(unittest.TestCase):
    
    def setUp(self):
        pass
    
    def test_rect(self):
        r = cl2shp.Rect( [-100,80], [-60,84] )
        self.assertTrue ( r.contains( ( -80,82) ) )
        self.assertFalse( r.contains( (-100,82) ) )
        self.assertFalse( r.contains( ( -80,79) ) )
        self.assertTrue ( r.contains( ( -90.5,83.4) ) )
    
    def test_connect_lines(self):
        island_way_id = 10663437
        
        chopper = cl2shp.CoastlineChopper(test_data_file, rectangle, threads=num_threads)
        
        island_len = len(chopper.coastlines[island_way_id].nodes())
        
        chopper.connect_lines()
        
        # check if islands are unchanged
        self.assertEqual( island_len, len(chopper.coastlines[island_way_id].nodes()) )
        
        ways_a = getWaysContainingNodeId(chopper.coastlines.viewvalues(), 2716005436)
        ways_b = getWaysContainingNodeId(chopper.coastlines.viewvalues(), 2146092568)

        # check that both lists contain the same ways
        self.assertEqual( ways_a[0], ways_b[0] )
        # check that the node is only contained in one coastline
        self.assertEqual( 1, len(ways_a) )
        self.assertEqual( 1, len(ways_b) )

    def test_chop_ways(self):
        chopper = cl2shp.CoastlineChopper(test_data_file, rectangle, threads=num_threads)
        chopper.connect_lines()
        chopper.chop_ways()
        
        # assert that we still have coordindates
        coordinates = chopper.coordinates.viewvalues()
        self.assertTrue( len(coordinates) > 0 )
        
        # assert that all nodes are contained in the rectangle
        boarder_node_ids = chopper.boarder_nodes.keys()
        for coord in coordinates:
            self.assertTrue( rectangle.contains(coord[1:3]) or (coord[0] in boarder_node_ids) )
        
    def test_close_coastlines(self):
        cl_util = cl2shp.CoastlineChopper(test_data_file, rectangle, threads=num_threads)
        
        cl_util.connect_lines()
        cl_util.chop_ways()
        cl_util.close_open_lines()
        
        for coastline in cl_util.coastlines.viewvalues():
            # check the coastline has been closed
            self.assertTrue( coastline.is_circle() )
            
            # check the direction of the coastline is anti-clockwise around land
            doubled_area = 0
            last_coord = cl_util.coordinates[coastline.last_node()][1:3]
            for nid in coastline.nodes():
                coord = cl_util.coordinates[nid][1:3]
                doubled_area += last_coord[0] * (coord[1] - last_coord[1])
                last_coord = coord
            
            self.assertTrue( doubled_area > 0 )
            
    
if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(cl2shpTester)
    unittest.TextTestRunner(verbosity=2).run(suite)
