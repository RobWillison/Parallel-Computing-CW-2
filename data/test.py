processors = [1, 2, 4, 8, 16, 32, 48, 64]

fout = open('efficiency9000','w')

first = True

with open("TimeThread9000", "r") as my_file:
    for lines in my_file:
        lines = lines.replace("\n", '');
        numbers = lines.split(' ')

        if first:
            first = False
            continue

        if int(numbers[0]) == 1:
            core1 = float(numbers[1])
            continue

        p = float(numbers[0])

        result = core1 / (p * float(numbers[1]))

        fout.write(str(int(p)) + ' ' + str(result) + '\n')

        print(p)

fout.close()
