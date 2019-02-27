def make_csv(infile, output, separator=',', add_string=False):
    """
    Takes an input file, strips it of all empty spaces and turns it into a comma separated file
    @type infile: file name
    @type output: file name
    @param infile: name of input file
    @param output: name of output file
    @type separator: String
    @param add_string: Field Separator for output file (defaul is comma)
    @type add_string: False/String
    @param add_string: Enter any string you might want to add to the beginning of each line, default is False
    @note: This function will happily overwrite any important file you give as output!
    """
    print 'Input File is: ' + infile
    print 'output File is: ' + output
    if add_string:
        print 'Adding' + add_string
    f = open(infile, 'r')
    out = open(output, 'w')
    for line in f.readlines():
        tmp = ''
        if line[0] == '#':
            pass
        else:
            if add_string:
                out.write(add_string)
            for x in line.split():
                tmp += x.strip() + separator
            out.write(tmp[:-1] + '\n')
    out.close()
    print "Done!"


