from imposm.parser import OSMParser
import codecs


output_file = 'data/sweden-coastline.osm'
boarders = ((12.0, 55.0), (14.0, 56.0))


# multidimensional implementation of line segments
class Line:
    def __init__(self, origin, destination):
        self.origin = origin
        self.dest   = destination
    
    # interpolate the coordinate 'p' for 'p[dim] = val'
    def interpolate(self, dim, val):
        result = []
        
        q = (val - self.origin[dim]) / (self.dest[dim] - self.origin[dim])
        
        for d in range( len(self.origin) ):
            result.append( val if d == dim else q * (self.dest[d] - self.origin[d]) + self.origin[d] )
        
        return tuple(result)
    
    def intersect(self, other):
        ds = self.vector()
        do = other.vector()
        
        fs  = float(self.origin[0] - other.origin[0]) / do[0] * (float(ds[1]) / do[1]) + float(other.origin[1] - self.origin[1]) / ds[1]
        fs /= (1 - float(ds[0] * do[1]) / (do[0] * ds[1]))
        
        fo = (fs * ds[0] + self.origin[0] - other.origin[0]) / do[0]
        
        return (fs, fo)
    
    
    def vector(self):
        d = []
        for dim in range( len(self.origin) ):
            d.append(self.dest[dim] - self.origin[dim])
        
        return tuple(d)
    

class Polygon:
    def __init__(self, points):
        self.points = points
    
    def number_of_edges(self):
        return len(self.points) - 1
    
    def contains(self, coord):
        # see: https://en.wikipedia.org/wiki/Point_in_polygon#Ray_casting_algorithm
        
        crossings = 0
        last_point = self.points[0]
        
        for point in self.points:
            if last_point[0] < coord[0] < point[0] or last_point[0] > coord[0] > point[0]:
                l = Line(last_point, point)
                pos = l.interpolate(0, coord[0])
                
                if pos[1] > coord[1]:
                    crossings += 1
        
        return (crossings % 2) == 1
    
    
    def intersect_with_lineseg(self, line):
        
        last_point = None
        
        for point in self.points:
            if last_point is None:
                last_point = point
                continue
            
            l = Line(last_point, point)
            f = l.intersect(line)
            
            if not (0 <= f[0] <= 1 and 0 <= f[1] <= 1):
                continue
            else:
                break
        
        x = last_point[0] + float(point[0] - last_point[0]) * f[0]
        y = last_point[1] + float(point[1] - last_point[1]) * f[1]
        
        return (x, y)
    
    

class Rect:
    def __init__(self, origin, other):
        self.origin = origin
        self.other = other
    
    def contains(self, coord):
        result = True
        
        for d in range( len(self.origin) ):
            result = result and ( self.origin[d] <= coord[d] <= self.other[d] or self.origin[d] >= coord[d] >= self.other[d])
        
        return result
    
    def next_coord_on_boarder(self, coord):
        
        if not self.contains( coord ):
            return None
        
        dims = len(self.origin)
        min_dist = abs(self.other[0] - self.origin[0])
        edge = -1
        res_coord = None
        
        for d in range( dims ):
            
            dist_orig = abs(coord[d] - self.origin[d])
            
            if dist_orig < min_dist:
                min_dist = dist_orig
                edge = dims * (d % 2) + d
            
            dist_other = abs(coord[d] - self.other[d])
            
            if dist_other < min_dist:
                min_dist = dist_other
                edge = dims * ((d + 1) % 2) + d
        
        if edge == -1:
            return None
        
        res_coord = list(coord)
        # get the dimension of the closest edge
        d = edge % dims
        # project the coordinate on the corresponding edge
        res_coord[d] = self.origin[d] if d == edge else self.other[d]
        
        return (edge, tuple(res_coord))
    
    
    def intersect_with_lineseg(self, line):
        
        result = None
        edge = -1
        
        dims = len(self.origin)
        
        for d in range( dims ):
            if   line.origin[d] < self.origin[d] < line.dest[d] or line.origin[d] > self.origin[d] > line.dest[d]:
                result = line.interpolate(d, self.origin[d])
                edge = dims * (d % 2) + d
            elif line.origin[d] < self.other[d]  < line.dest[d] or line.origin[d] > self.other[d]  > line.dest[d]:
                result = line.interpolate(d, self.other[d])
                edge = dims * ((d + 1) % 2) + d
            
            # update line.origin or line.dest with result for the case that line intersects more than one of the boarders of the rectangle, outside of it
            if result is not None:
                if not self.contains(line.origin):
                    line = Line(result, line.dest)
                else:
                    line = Line(line.origin, result)
        
        return (edge, result) if edge != -1 else None
    
    def number_of_edges(self):
        return 4
    
    
    def corner(self, index):
        i = index % 4
        
        if i == 0:
            return tuple(self.origin)
        elif i == 1:
            return (self.origin[0], self.other[1])
        elif i == 2:
            return tuple(self.other)
        else:
            return (self.other[0], self.origin[1])
        


class OSMWay(object):
    def __init__(self, imposm_way=None, osm_id=None):
        if imposm_way is None:
            self.tags = {}
            self.refs = []
        else:
            self.id = imposm_way[0]
            self.tags = imposm_way[1]
            self.refs = imposm_way[2]
        
        if osm_id is not None:
            self.id = osm_id
        
    
    def get_id(self):
        return self.id
    
    def __len__(self):
        return len(self.refs) - 1
    
    def append(self, other_way):
        try:
            self.refs.extend(
                other_way.refs if (self.refs[-1] != other_way.refs[0]) else other_way.refs[1:]
            )
        except IndexError:
            # the first or the second way may be empty
            self.refs.extend(other_way.refs)
    
    def append_node(self, ref):
        self.refs.append(ref)
    
    def first_node(self):
        return self.refs[0]
    
    def last_node(self):
        return self.refs[-1]
    
    def nodes(self):
        return self.refs
    
    def is_circle(self):
        return self.refs[0] == self.refs[-1]
    
    def insert_node(self, index, ref):
        self.refs.insert(index, ref)
    
    # creates a new way using all nodes between start_node and end_node (start_node and end_node inclusive)
    def extract_part(self, start_node, end_node, include_last=False, osm_id=0):
        start_pos = self.refs.index(start_node)
        end_pos = self.refs.index(end_node, start_pos + 1) + (1 if include_last else 0)
        
        return OSMWay( (osm_id, self.tags, self.refs[start_pos:end_pos]) )
    
    def remove_ind(self, i):
        del self.refs[i]
    
    def remove_node(self, ref):
        self.refs.remove(ref)
    
    def to_shape(self, coord_dict):
        shape = [ list(coord_dict[ref][1:3]) for ref in self.refs ]
        shape.reverse()
        
        return shape



# simple class that handles the parsed OSM data.
class CoastlineChopper(object):
    
    def __init__(self, input_file, boarder_rect):
        self.input_file = input_file
        self.boarder_rect = boarder_rect
        
        self.coordinates = {}
        self.coastlines  = {}
        self.continent_ids = set()
        self.boarder_nodes = {}
        self.auxiliary_items = 0L
        
        self.load_data()
        
    
    def load_data(self):
        
        print('loading data from OSM-file "%s"...' % self.input_file)
        
        self.coordinates = {}
        self.coastlines  = {}
        
        parser = OSMParser( concurrency=4,
                            coords_callback = self._load_coords_callback,
                            ways_callback   = self._load_coastlines_callback )
        parser.parse(self.input_file)
        
        return
    
    def _load_coords_callback(self, coords):
        # build a hash-table, hashing the osmids on the coordinates
        for c in coords:
            self.coordinates[ c[0] ] = c
        return
    
    def _load_coastlines_callback(self, ways):
        # callback method for ways
        for osmid, tags, refs in ways:
            if 'natural' in tags and tags['natural'] == 'coastline':
                self.coastlines[osmid] = OSMWay( imposm_way=(osmid, tags, refs) )
        return
    
    
    def get_node_ids(self):
        print('collecting node ids...')
        
        node_ids = []
        
        for way in self.coastlines.viewvalues():
            node_ids.extend(way.nodes())
        return set(node_ids)
    
    
    def connect_lines(self):
        print('connecting lines...')
        
        first_nodes = {}
        
        for osmway in self.coastlines.viewvalues():
            first_nodes[osmway.first_node()] = osmway
        
        # blacklist containing the ids of the ways, that have already been appended to others
        blacklist = set()
        
        for osmway in self.coastlines.viewvalues():
            
            if osmway.get_id() in blacklist:
                continue
            
            try:
                while True:
                    other_way = first_nodes[osmway.last_node()]
                    
                    if other_way is osmway:
                        # reached a circle
                        break
                    
                    osmway.append( other_way )
                    
                    # blacklist the way that has just been appended to the other
                    blacklist.add(other_way.get_id())
                
            except KeyError:
                # line ends and has no sucessor
                pass
        
        # clean up
        print('removing %d lines not needed any more...' % len(blacklist))
        for way_id in blacklist:
            del self.coastlines[way_id]
        
        return
    
    def _get_next_id(self):
        self.auxiliary_items += 1
        return (-1) * self.auxiliary_items
    
    def _create_node(self, coords):
        nid = self._get_next_id()
        node = (nid, coords[0], coords[1])
        
        self.coordinates[nid] = node
        return nid
    
    
    def line_stats(self):
        circles = 0
        open_lines = 0
        open_lines_lengths = []
        
        for osmway in self.coastlines.viewvalues():
            if osmway.is_circle():
                circles += 1
            else:
                open_lines += 1
                open_lines_lengths.append( len(osmway) )
        
        return (circles, open_lines, open_lines_lengths)
    
    
    def _add_boarder_node(self, direction, coord1, coord2, osmway):
        # create the id for the node on the boarder
        bn_id = self._get_next_id()
        # and append it to the new way
        if direction == 'in':
            osmway.insert_node(0, bn_id)
        else:
            osmway.append_node(bn_id)
        
        if coord1 is not None and coord2 is not None:
            # calculate intersection (Line coord1, coord2 with the boarder_rect)
            # 1. the edge that is intersected by the line
            # 2. the coordinate of the intersection
            (edge_id, coord) = self.boarder_rect.intersect_with_lineseg(
                # coord[0] is the id, so we start at 1
                Line(coord1, coord2)
            )
        else:
            # calculate the next coordinate to coord1 on the rectangle
            (edge_id, coord) = self.boarder_rect.next_coord_on_boarder(
                coord1 if coord1 is not None else coord2
            )
        
        bn = (direction, bn_id, osmway)
        
        self.boarder_nodes[edge_id].append( bn )
        
        return (bn_id, coord[0], coord[1])
    
    
    def chop_ways(self):
        self.boarder_nodes = self.boarder_rect.number_of_edges() * [ [] ]
        
        # rebuild the continent ids
        self.continent_ids = set()
        
        filtered_coords = {}
        chopped_coastlines = {}
        
        print('building boundaries')
        
        for osmway in self.coastlines.viewvalues():
            
            beginning_inside = False
            first_seg_id = None
            
            all_nodes_inside = True
            entrance_ref = None
            cur_id = None
            
            last_coord_inside = None
            last_coord = None
            coord_inside = None
            coord = None
            
            is_circle = osmway.is_circle()
            
            for ref in osmway.nodes():
                
                last_coord_inside = coord_inside
                last_coord = coord
                
                try:
                    coord = self.coordinates[ref][1:]
                except KeyError:
                    print('\tWARNING: undefined reference %d in way %d' % (ref, osmway.get_id()))
                    continue
                
                coord_inside = self.boarder_rect.contains(coord)
                
                if coord_inside:
                    # add coordinate to our new coordinate list
                    filtered_coords[ref] = (ref, coord[0], coord[1])
                    
                    if last_coord_inside:
                        continue
                    
                    entrance_ref = ref
                    
                    # create the new OSMWay for this piece of coastline
                    cur_id = self._get_next_id()
                    chopped_coastlines[cur_id] = OSMWay(osm_id=cur_id)
                    
                    if last_coord_inside is None:
                        # this is the first coordinate of the current way
                        # and it is inside
                        # this makes sense for islands, for others this is bad
                        beginning_inside = True
                        first_seg_id = cur_id
                        
                        if is_circle:
                            # do not insert a boarder node!
                            # this line is a circle, and the line will be extended
                            # in the front later
                            continue
                        
                        print('WARNING: coastline with id %d is beginning inside!' % osmway.get_id())
                        
                    else:
                        print('\tline %d going inside' % osmway.get_id())
                
                else:
                    # coordinate is outside of the area we want to extract
                    all_nodes_inside = False
                    
                    if not last_coord_inside:
                        continue
                    
                    # coastline leaves our area of interest
                    print('\tline %d going outside' % osmway.get_id())
                    
                    # this is a continental line
                    self.continent_ids.add(cur_id)
                    
                    # create a new way from the last piece of this osmway that has been inside
                    chopped_coastlines[cur_id].append(
                        osmway.extract_part(entrance_ref, ref)
                    )
                
                # calculate the boarder node
                bn = self._add_boarder_node(
                    'in' if coord_inside else 'out',
                    last_coord,
                    coord,
                    chopped_coastlines[cur_id]
                )
                filtered_coords[bn[0]] = bn
            
            # finish the line, if the last node is still inside
            if coord_inside:
                chopped_coastlines[cur_id].append(
                    osmway.extract_part(entrance_ref, ref, include_last=True)
                )
                
                if not is_circle:
                    print('WARNING: coastline with id %d is ending inside!' % osmway.get_id())
                    
                    # get the id for the boarder node
                    bn = self._add_boarder_node('out', coord, None, chopped_coastlines[cur_id])
                    filtered_coords[bn[0]] = bn
                    
                elif not all_nodes_inside:
                    # merge the last and the first piece
                    chopped_coastlines[cur_id].append(chopped_coastlines[first_seg_id])
                    
                    del chopped_coastlines[first_seg_id]
                    
                    # fix the continent ids
                    self.continent_ids.remove(first_seg_id)
                    self.continent_ids.add(cur_id)
            
        
        print('filtered out %d of %d nodes and %d of %d coastlines' %
              (len(filtered_coords),
               len(self.coordinates),
               len(chopped_coastlines),
               len(self.coastlines)) )
        
        # delete old collection of coordinates
        del self.coordinates
        # and replace it by the new one
        self.coordinates = filtered_coords
        
        # tidy up the coastlines
        del self.coastlines
        self.coastlines = chopped_coastlines
        
        return self.boarder_nodes
    
    
    def number_of_coastlines(self):
        return len(self.coastlines)
    
    
    def close_open_lines(self):
        
        sorted_boarder_nodes = []
        for edge in range( self.boarder_rect.number_of_edges() ):
            # as long as we are using rectangles this works
            # has to be fixed for polygons
            direction = -1 if edge >= 2 else 1
            f = lambda bn: direction * self.coordinates[ bn[1] ][(edge % 2) + 1]
            
            sorted_boarder_nodes.append( sorted( self.boarder_nodes[edge], key=f ) )
        
        blacklist = set()
        
        first_way = None
        inside = None
        corner_nodes = []
        
        for edge in range( len(sorted_boarder_nodes) ):
            corner_nodes.append( self._create_node(self.boarder_rect.corner(edge)) )
            
            for bn in sorted_boarder_nodes[edge]:
                if bn[0] == 'out':
                    while corner_nodes:
                        # append the corners of the rectangle/polygon that have passed
                        bn[2].append_node( corner_nodes.pop() )
                    
                    if in_node is not None:
                        # connect this coastline with the next one
                        if in_node[2] is bn[2]:
                            # circle
                            bn[2].append_node( bn[2].first_node() )
                        else:
                            # other coastline
                            bn[2].append( in_node[2] )
                            blacklist.add( in_node[2].get_id() )
                            
                            if first_way is in_node[2]:
                                first_way = bn[2]
                    else:
                        first_way = bn[2]
                    
                    in_node = None
                    
                else:
                    in_node = bn
                    corner_nodes = []
        
        if first_way is not None:
            while corner_nodes:
                first_way.append_node( corner_nodes.pop() )
            
            first_way.append_node( first_way.first_node() )
        
        for way_id in blacklist:
            del self.coastlines[way_id]
        
        return
    
    def check_coordinates(self):
        for osmway in self.coastlines.viewvalues():
            for ref in osmway.nodes():
                try:
                    self.coordinates[ref]
                except KeyError:
                    print('ERROR: node %d missing!' % ref)
    
    def print_local_nodes(self):
        for node_id in filter(lambda i: i<0, self.coordinates.keys()):
            print(self.coordinates[node_id])
    
    
    def write_shapefile(self, output_file):
        
        # import shapefile module
        import shapefile
        
        w = shapefile.Writer(shapefile.POLYGON)
        w.autoBalance = 1
        
        w.field('OSMID')
        w.field('type', 'C', '15')
        
        for coastline in self.coastlines.viewvalues():
            print('%d: %d nodes' % (coastline.get_id(), len(coastline)))
            w.record( str(coastline.get_id()), 'coastline' )
            w.poly( shapeType=5, parts=[coastline.to_shape(self.coordinates)] )
        
        w.save( output_file )
        
        return
    

if __name__ == '__main__':
    
    import sys
    
    input_file = sys.argv[1]
    boarder_rect = Rect( boarders[0], boarders[1] )
    
    # instantiate counter and parser and start parsing
    cl_util = CoastlineChopper( input_file, boarder_rect )
    
    print( '%d coastlines found\n' % cl_util.number_of_coastlines() )
    cl_util.connect_lines()
    
    print( 'filter ways and coordinates in the region...' )
    cl_util.chop_ways()
    
    print( '%d coastlines after merging' % cl_util.number_of_coastlines() )
    print('islands: %d\ncontinents: %d\n' % cl_util.line_stats()[0:2])
    
    cl_util.close_open_lines()
    
    print( 'writing shapefile \"%s\"' % 'sweden-coastlines.shp' )
    cl_util.write_shapefile('sweden-coastlines.shp')
    
