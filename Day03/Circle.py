import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--radius', help='Enter radius', required=True, type=int)

args = parser.parse_args()

radius = args.radius

area = 3.14*radius*radius
circumference = 2*3.14*radius
print(" Area of circle is:", area, '\n', "Circumference of circle is:", circumference)