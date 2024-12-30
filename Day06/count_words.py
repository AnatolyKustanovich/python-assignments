def count_words(text):
    words = text.split()
    return len(words)

def count_number_of_each_word(text):
    words = text.split(",")
    word_count = {}
    for word in words:
        if word in word_count:
            word_count[word] +=1
       #word_count[word] = word_count[word] + 1
        else:
            word_count[word] = 1
    return word_count

print(count_words("name,23,value,42.1"))
res = print(count_number_of_each_word("name,23,value,42.1,name,23"))
print(res)
for key in res.keys():
    print(key, res[key])
