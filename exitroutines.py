def getout(numberofwarnings, exitcode):
    if exitcode < 0:
        print('Abnormal interruption')
        print('WARNINGS: ', numberofwarnings)
        raise SystemExit
    elif exitcode > 0:
        print('Success!')
        print('WARNINGS: ', numberofwarnings)
