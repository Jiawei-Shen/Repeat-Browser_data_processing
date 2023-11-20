/* paraFetch - fetch things from remote URLs in parallel. */

#ifndef PARAFETCH_H
#define PARAFETCH_H

boolean parallelFetch(char *url, char *outPath, int numConnections, int numRetries, boolean newer, boolean progress);
/* Open multiple parallel connections to URL to speed downloading */

struct parallelConn
/* struct to information on a parallel connection */
    {
    struct parallelConn *next;  /* next connection */
    int sd;                     /* socket descriptor */
    off_t rangeStart;           /* where does the range start */
    off_t partSize;             /* range size */
    off_t received;             /* bytes received */
    };

boolean paraFetchReadStatus(char *origPath, 
    struct parallelConn **pPcList, char **pUrl, off_t *pFileSize, 
    char **pDateString, off_t *pTotalDownloaded);
/* Read a status file - which is just origPath plus .paraFetchStatus.  This is updated during 
 * transit by parallelFetch. Returns FALSE if status file not there - possibly because
 * transfer is finished.  Any of the return parameters (pThis and pThat) may be NULL */

#endif /* PARAFETCH_H */
