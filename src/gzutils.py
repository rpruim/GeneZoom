#!/usr/bin/env python

#===============================================================================
# Copyright (C) 2010 Ryan Welch, Randall Pruim
# 
# LocusZoom is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# LocusZoom is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#===============================================================================

import sys
import os
import decimal
import fnmatch
from which import which_files

def check(expression, msg=None):
	if msg == None:
		pass
	else:
		print "\n", msg, "... ",
	try:
		assert(expression)
		print "Pass ",
		return(True)
	except:
		print "Fail ",
		return(False)

def find_relative(file):
  full_path = None;
  
  # Find the m2zfast root, using the script's location. 
  start_loc = os.path.realpath(sys.argv[0]);
  script_dir = None;
  if os.path.isdir(start_loc):
    script_dir = start_loc;
  else:
    script_dir = os.path.dirname(start_loc);
  root_dir = os.path.join(script_dir,"../");
  
  # If the file to find has a path, it means it is a path relative
  # to the m2zfast root. We need to attach that path to the root. 
  (file_path,file_name) = os.path.split(file);
  if file_path != "":
    root_dir = os.path.abspath(os.path.join(root_dir,file_path));
  
  if file_name == "" or file_name == None:
    if os.path.exists(root_dir):
      full_path = root_dir;
  else:
    temp_path = os.path.join(root_dir,file_name);
    if os.path.exists(temp_path):
      full_path = temp_path;
  
  return full_path;

def locate(pattern, root=os.curdir, followlinks=True):
  '''Locate all files matching supplied filename pattern in and below
  supplied root directory.'''
  for path, dirs, files in os.walk(os.path.abspath(root),followlinks=followlinks):
    for dir in fnmatch.filter(dirs,pattern):
      yield os.path.join(path,dir);
    for filename in fnmatch.filter(files, pattern):
      yield os.path.join(path,filename)

# Tries to find a file in the following order:
# 1) the given file (duh)
# 2) relative to m2zfast's root directory
# 3) on the user's path
def find_systematic(file):
  if os.path.isfile(file):
    return os.path.abspath(file);

  relative = find_relative(file);
  if relative:
    return relative;
  
  walk = [i for i in which_files(file)];
  if len(walk) > 0:
    return walk[0];
    
  return None;

def is_number(x):
  try:
    decimal.Decimal(str(x));
  except:
    return False;

  return True;

# Terminates program with only a message.
def die(msg):
  print >> sys.stderr, msg;
  sys.exit(1);

# Continue program after warning message.
def warn(msg, e=None):
  print >> sys.stderr, msg;
  if e:
	  print >> sys.stderr, e
  print >> sys.stderr, "\n"

# Tests to see if interval i1 is completely contained within interval i2.
# i1 should be a list or tuple [start,end], same with i2
def interval_contained(i1,i2):
  if i1[0] >= i2[0] and i1[1] <= i2[1]:
    return True;
  else:
    return False;

def chrom2chr(c):
  c = str(c);
  if c in ('X','chrX','chromX'):
    return 23;
  elif c in ('Y','chrY','chromY'):
    return 24;
  elif c == 'mito':
    return 25;
  elif c == 'XY':
    return 26;
  else:
    try:
      c = int(c);
    except:
      c = None;

    return c;

def chr2chrom(c):
  c = int(c);
  if c < 23 and c > 0:
    return str(c);
  elif c == 23:
    return 'X';
  elif c == 24:
    return 'Y';
  elif c == 25:
    return 'mito';
  elif c == 26:
    return 'XY';
  else:
    return None;
  
