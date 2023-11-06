import h5py
import sys

# FILE ~ *.hdf5
FILE = sys.argv[1]

def main(FILE):
  f = h5py.File(FILE,'a')
  #['atomic_masses', 'cartesian_coordinates', 'elements', ~]
  config = f['BulkConfiguration_0']['AtomicConfiguration']
  ele  = config['elements']['atomic_numbers']['data']
  #mass = config['atomic_masses']['array']['data']
  #pos  = config['cartesian_coordinates']['array']['data']

  count = count_ions(ele)

###
def count_ions(ELEMENTS):
  count = {}
  for i in range(len(ELEMENTS)):
    if ELEMENTS[i] not in count.keys():
      count[ELEMENTS[i]] = 1
    else:
      count[ELEMENTS[i]] += 1
  for idx, item in enumerate(count):
    print(item, count[item])
  return count
  
###
if __name__ == "__main__":
    main(FILE)
