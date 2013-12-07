#!/usr/bin/python2
# -*- coding: utf-8 -*-

# Copyright © 2013 Raphael Dümig <duemig@in.tum.de>
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# 
# 
# Dieses Programm ist Freie Software: Sie können es unter den Bedingungen
# der GNU General Public License, wie von der Free Software Foundation,
# Version 3 der Lizenz oder (nach Ihrer Wahl) jeder neueren
# veröffentlichten Version, weiterverbreiten und/oder modifizieren.
# 
# Dieses Programm wird in der Hoffnung, dass es nützlich sein wird, aber
# OHNE JEDE GEWÄHRLEISTUNG, bereitgestellt; sogar ohne die implizite
# Gewährleistung der MARKTFÄHIGKEIT oder EIGNUNG FÜR EINEN BESTIMMTEN ZWECK.
# Siehe die GNU General Public License für weitere Details.
# 
# Sie sollten eine Kopie der GNU General Public License zusammen mit diesem
# Programm erhalten haben. Wenn nicht, siehe <http://www.gnu.org/licenses/>. 


# python2 is needed for the dependency imposm
from imposm.parser import OSMParser
import codecs

# if you want no color in the shell by default, set this to 'True'
color_disabled = False

# message highlighting
shellColor = { 'black':  30,
               'red':    31,
               'green':  32,
               'yellow': 33,
               'blue':   34,
               'purple': 35,
               'cyan':   36,
               'white':  37,
               'off':     0
             }
shellBold  = { False: 0,
               True:  1
             }


def styled(text, color='off', bold=False):
    if color_disabled:
        return text
    elif color is not None:
        return ('\033[%d;%dm' % (shellBold[bold], shellColor[color])) + text + ('\033[%dm' % shellColor['off'])



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
    

# not used so far: intended to replace Rect at some time
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
    
    def __str__(self):
        return '\t%f\n%f\t%f\n\t%f\n' % (self.other[1], self.origin[0], self.other[0], self.origin[1])
    
    def contains(self, coord):
        result = True
        
        for d in range( len(self.origin) ):
            result = result and ( self.origin[d] < coord[d] < self.other[d] or self.origin[d] > coord[d] > self.other[d])
        
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


# the main class of this script
class CoastlineChopper(object):
    
    def __init__(self, input_file, boarder_rect, threads=1, node_filter=True):
        self.input_file = input_file
        self.boarder_rect = boarder_rect
        
        self.coordinates = {}
        self.coastlines  = {}
        self.continent_ids = set()
        self.boarder_nodes = {}
        self.auxiliary_items = 0L
        
        self.threads = threads
        
        self.load_data(node_filter)
    
    
    def load_data(self, node_filter=True):
        
        self.coordinates = {}
        self.coastlines  = {}
        
        if not node_filter:
            # load all data at once
            # problem: we have to save a hashtable with ALL coordinates
            parser = OSMParser( concurrency     = self.threads,
                                coords_callback = self._load_coords_callback,
                                ways_callback   = self._load_coastlines_callback )
            parser.parse(self.input_file)
            
        else:
            # first load the coastlines
            way_parser = OSMParser( concurrency = self.threads,
                                    ways_callback = self._load_coastlines_callback )
            way_parser.parse(self.input_file)
            
            # now load only the coordinates used in the coastlines
            node_ids = self.get_node_ids()
            
            print('\tloading coordinates...')
            # this saves an amount of memory, but is much slower
            # as we have to look up each coordinate
            coord_parser = OSMParser( concurrency = self.threads,
                                      coords_callback = lambda coords: self._load_coords_callback(coords, filter_set=node_ids) )
            coord_parser.parse(self.input_file)
        
        return
    
    def _load_coords_callback(self, coords, filter_set=None):
        # build a hash-table, hashing the OSM-ids on the coordinates
        
        coord_list = coords if filter_set is None else filter(lambda c: c[0] in filter_set, coords)
        
        for c in coord_list:
            self.coordinates[ c[0] ] = c
        
        return
    
    def _load_coastlines_callback(self, ways):
        # callback method for ways
        for osmid, tags, refs in ways:
            if 'natural' in tags and tags['natural'] == 'coastline':
                self.coastlines[osmid] = OSMWay( imposm_way=(osmid, tags, refs) )
        return
    
    
    def get_node_ids(self):
        print('\tcollecting node ids...')
        
        node_ids = set()
        
        for way in self.coastlines.viewvalues():
            for node in way.nodes():
                node_ids.add(node)
        
        return node_ids
    
    
    def connect_lines(self):
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
                # => coastline of a continent
                pass
        
        # clean up
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
                
                if self.boarder_rect.contains(self.coordinates[osmway.first_node()][1:3]):
                    print('WARNING: open line %d beginning inside!' % osmway.get_id())
                if self.boarder_rect.contains(self.coordinates[osmway.last_node()][1:3]):
                    print('WARNING: open line %d ending inside!' % osmway.get_id())
        
        return (circles, open_lines)
    
    
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
        
        bn = (bn_id, edge_id, direction, osmway)
        
        self.boarder_nodes[bn_id] = bn
        
        return (bn_id, coord[0], coord[1])
    
    
    def chop_ways(self):
        self.boarder_nodes = {}
        
        # rebuild the continent ids
        self.continent_ids = set()
        
        filtered_coords = {}
        chopped_coastlines = {}
        
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
                    coord = self.coordinates[ref][1:3]
                except KeyError:
                    print('\tWARNING: undefined node reference %d in way %d' % (ref, osmway.get_id()))
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
                        # and it is already inside of our boundary rectangle
                        # 
                        # this makes sense for islands, as the coastline is continued
                        # at the other end of the line, but for others this is bad,
                        # as we have to find some way to connect it to the boundary
                        # of our boundary rectangle
                        beginning_inside = True
                        first_seg_id = cur_id
                        
                        if is_circle:
                            # do not insert a boarder node!
                            # this line is a circle, and the line will be extended
                            # in the front later
                            continue
                        else:
                            print('\tWARNING: coastline with id %d is beginning inside!' % osmway.get_id())
                
                else:
                    # coordinate is outside of the area we want to extract
                    all_nodes_inside = False
                    
                    if not last_coord_inside:
                        continue
                    
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
                
            # all nodes of this way have been parsed by now
            # 
            # finish the line, if the last node is still inside
            if coord_inside:
                chopped_coastlines[cur_id].append(
                    osmway.extract_part(entrance_ref, ref, include_last=True)
                )
                
                if not is_circle:
                    print('\tWARNING: coastline with id %d is ending inside!' % osmway.get_id())
                    
                    # get the id for the boarder node
                    bn = self._add_boarder_node('out', coord, None, chopped_coastlines[cur_id])
                    filtered_coords[bn[0]] = bn
                    
                elif not all_nodes_inside:
                    # circle, but not all nodes inside:
                    # we once had a circle, but by now it has been separated to pieces
                    # 
                    # so we append the first to the last piece
                    first_cl = chopped_coastlines[first_seg_id]
                    
                    chopped_coastlines[cur_id].append(first_cl)
                    
                    # find the boarder node of the first segment, and link it to the current/last segment
                    bn = list(self.boarder_nodes[first_cl.last_node()])
                    bn[-1] = chopped_coastlines[cur_id]
                    self.boarder_nodes[first_cl.last_node()] = bn
                    
                    del chopped_coastlines[first_seg_id]
                    
                    # fix the continent ids
                    self.continent_ids.remove(first_seg_id)
                    self.continent_ids.add(cur_id)
            
        
        # delete old collection of coordinates
        del self.coordinates
        # and replace it by the new one
        self.coordinates = filtered_coords
        
        # tidy up the coastlines
        del self.coastlines
        self.coastlines = chopped_coastlines
        
        return self.boarder_nodes
    
    
    def close_open_lines(self):
        
        # sort the boarder nodes: first by edges, then by coordinates
        
        # initialize sorted list
        sorted_boarder_nodes = []
        for edge_id in range(self.boarder_rect.number_of_edges() ):
            sorted_boarder_nodes.append( [] )
        
        # sort by edges
        for bn in self.boarder_nodes.viewvalues():
            sorted_boarder_nodes[ bn[1] ].append(bn)
        
        for edge_id in range( self.boarder_rect.number_of_edges() ):
            # sorting by coordinates
            
            # as long as we are using rectangles this works
            # but this has to be fixed for polygons
            direction = -1 if edge_id >= 2 else 1
            f = lambda bn: direction * self.coordinates[ bn[0] ][(edge_id + 1) % 2 + 1]
            
            sorted_boarder_nodes[edge_id] = sorted( sorted_boarder_nodes[edge_id], key=f )
        
        
        blacklist = set()
        aliases = {}
        first_way = None
        in_node = None
        inside = None
        corner_nodes = []
        
        
        for edge in range( len(sorted_boarder_nodes) ):
            corner_nodes.append( self._create_node(self.boarder_rect.corner(edge)) )
            
            for bn in sorted_boarder_nodes[edge]:
                if bn[2] == 'out':
                    cl = bn[-1]
                    
                    # if the coastline of bn has already been replaced by another one, and this again
                    # by another one and so on, load this instead of the original coastline
                    # else simply use the coastline of bn
                    try:
                        while True:
                            cl = aliases[cl]
                    except KeyError:
                        pass
                    
                    if corner_nodes:
                        print(styled('\tinserting %d corner nodes in way %d', color='yellow') % (len(corner_nodes), cl.get_id()))
                    
                    while corner_nodes:
                        # append the corners of the rectangle/polygon that have passed
                        cl.append_node( corner_nodes.pop() )
                    
                    if in_node is not None:
                        # connect this coastline with the next one
                        if in_node[-1] is cl:
                            # reached a circle
                            cl.append_node( cl.first_node() )
                            
                            print(styled('\tline %d closed', color='green') % cl.get_id())
                            
                        else:
                            # other coastline
                            other_cl = in_node[-1]
                            cl.append( other_cl )
                            blacklist.add( other_cl.get_id() )
                            aliases[other_cl] = cl
                            
                            if first_way is other_cl:
                                first_way = cl
                            
                            print(styled('\tlandmass out of boundary: appending %d to line %d', color='yellow') % (other_cl.get_id(), cl.get_id()))
                        
                        # remove the 'in_node'
                        in_node = None
                    else:
                        if first_way is None:
                            first_way = cl
                            print(styled('\tstarting position in landmass', color='yellow'))
                        else:
                            print(styled('ERROR: bad coastline (internal id=%d)', color='red') % cl.get_id())
                    
                    in_node = None
                    
                else:
                    in_node = bn
                    corner_nodes = []
        
        # if there is an open way left at the start and the end, connect the pieces
        if first_way is not None:
            while corner_nodes:
                first_way.append_node( corner_nodes.pop() )
            
            first_way.append_node( first_way.first_node() )
            print(styled('\tline %d closed', color='green') % first_way.get_id())
        
        # tidy up
        for way_id in blacklist:
            del self.coastlines[way_id]
        
        return
    
    
    def check_coordinates(self):
        for osmway in self.coastlines.viewvalues():
            for ref in osmway.nodes():
                try:
                    self.coordinates[ref]
                except KeyError:
                    print(styled('ERROR: node %d missing!', color='red') % ref)
    
    
    def write_shapefile(self, output_file):
        
        # import shapefile module
        import shapefile
        
        w = shapefile.Writer(shapefile.POLYGON)
        
        for coastline in self.coastlines.viewvalues():
            try:
                w.poly( shapeType=5, parts=[coastline.to_shape(self.coordinates)] )
            except KeyError:
                print('ERROR: undefined node in way %d: skipping!' % coastline.get_id())
        
        w.save( output_file )
        
        # generate the projection file
        basename = output_file[:-4] if output_file.endswith('.shp') else output_file
        with open(basename + '.prj', 'w') as prj:
            prj.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]')
        
        return
    

if __name__ == '__main__':
    
    # the command line
    import argparse
    
    parser = argparse.ArgumentParser(description='create a landmass shapefile from an OSM-extract')
    parser.add_argument('input_file', help='OSM-extract (XML/PBF) containing the coastlines')
    parser.add_argument('output_file', help='the destination for the shapefile')
    
    parser.add_argument('--bb', nargs=4, type=float, metavar=('NORTH', 'EAST', 'SOUTH', 'WEST'), default=[90.0,180.0,-90.0,-180.0], help='bounding box for the coastline data')
    parser.add_argument('--threads', type=int, default=2, help='maximum number of threads to use for loading the data (default: 2)')
    parser.add_argument('--no-node-filter', action='store_true', help='do not filter the nodes to those used by coastlines: All nodes will be loaded. Use this if you have already filtered the input data with for example osmosis. Results in faster parsing of the input data as only one pass through the data is needed.')
    parser.add_argument('--no-color', action='store_true', help='disable colorized output')
    
    args = parser.parse_args()
    
    boarder_rect = Rect( args.bb[3:1:-1], args.bb[1::-1] )
    if args.no_color:
        color_disabled = True
    
    prompt = styled('==> ', color='green', bold=True)
    print( prompt + styled('loading data from OSM-file "%s" (%d threads)', bold=True, color='white') % (args.input_file, args.threads))
    
    # instantiate counter and parser and start parsing
    cl_util = CoastlineChopper( args.input_file,
                                boarder_rect,
                                threads=args.threads,
                                node_filter=not args.no_node_filter )
    
    print( prompt + styled('connecting pieces of coastlines...', bold=True, color='white') )
    cl_util.connect_lines()
    print( '\tclosed coastlines: %d\topen coastlines: %d' % cl_util.line_stats() )
    
    print( prompt + styled('filtering ways and coordinates in the specified region', bold=True, color='white') )
    cl_util.chop_ways()
    print( '\tclosed coastlines: %d\topen coastlines: %d' % cl_util.line_stats() )
    
    print( prompt + styled('closing remaining coastlines...', bold=True, color='white') )
    cl_util.close_open_lines()
    print( '\tclosed coastlines: %d\topen coastlines: %d' % cl_util.line_stats() )
    
    print( prompt + styled('writing shapefile \"%s\"', bold=True, color='white') % args.output_file )
    cl_util.write_shapefile(args.output_file)
    
