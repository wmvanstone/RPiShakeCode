from datetime import datetime

def date2nthDay(date_string, format='%Y%m%d'):
    date_object = datetime.strptime(date_string, format)
    newYearDay = datetime(year=date_object.year, month = 1, day = 1)
    return(date_object - newYearDay).days + 1

def pad(num):
    out = str(num)
    while len(out) < 3:
        out = "0" + out
    return out

def bubblesort(ind, dist, col):
    for i in range(len(ind)-1,-1,-1):
        for j in range(i):
            if dist[ind[j]][col]>dist[ind[j+1]][col]:
                temp = ind[j]
                ind[j] = ind[j+1]
                ind[j+1] = temp
    return ind

if __name__ == "__main__":
    print("These utilities are designed to be run as a library")