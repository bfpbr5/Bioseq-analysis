import numpy as np

def readnpy(filename, outname):
    loadData=np.load(filename)  #加载文件
    doc = open(outname, 'a')  #打开一个存储文件，并依次写入
    # print(loadData, file=doc)  #将打印内容写入文件中
    print("----type----", file=doc)
    print(type(loadData), file=doc)
    print("----shape----", file=doc)
    print(loadData.shape, file=doc)
    print("----data----", file=doc)
    print(loadData, file=doc)
    doc.close()


def main():
    filename = 'genomic_benchmarks/my_datasets/pred_results.npy'
    outname = 'loaddata.txt'
    readnpy(filename,outname)


if __name__ == '__main__':
    main()