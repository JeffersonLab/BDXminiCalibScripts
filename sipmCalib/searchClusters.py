#This function, given the x and y positions of found peaks, determine clusters of nearby peaks and returns the SORTED array of peaks. The minimum distance to merge is thr.
def searchClusters(n,x,y,thr):
    #First, sort by peak height 
    xx=[]
    yy=[]
    clusters=[]
    for iii in range(n):
        xx.append(x[iii])
        yy.append(y[iii])
 
    yy,xx=(list(t) for t in zip(*sorted(zip(yy,xx),reverse=True)))

    #Now yy and xx are orderder descending wrt yy
    for iii in range(len(xx)):
        flag = False
        m_x1=xx[iii]
        m_y1=yy[iii]
        #Loop on all other hits
        for jjj in range(len(clusters)):
            m_x2=clusters[jjj]
            if (abs(m_x1-m_x2)<thr):                                            
                flag=True
                break
        
        if (flag == False):
            clusters.append(m_x1)
    
    clusters.sort()
    return clusters
