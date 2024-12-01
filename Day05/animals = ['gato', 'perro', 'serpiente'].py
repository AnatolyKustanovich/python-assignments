animals = ['gato', 'perro', 'serpiente']
#for animal in animals:
#    print(animal)

#animal = 'cat'
#while animal not in animals:
#    animal = input('type in the the name of the animal')

while True:
    animal = input('type in the the name of the animal')
    if animal in animals:
        break


text = "128746817 481"

count = 0
for char in text:
    if char.isdigit():
        count +=1
print(count)