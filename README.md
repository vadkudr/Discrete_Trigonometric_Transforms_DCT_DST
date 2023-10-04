**Discrete Trigonometric Transform**\
Fast methods for cosine and sine transforms computation through one another

| [Markus Puschel and Jose Moura The Algebraic Approach to the Discrete Cosine and Sine Transforms and their Fast Algorithms SIAM Journal of Computing 2003, Vol. 32, No. 5, pp. 1280-1316](http://www.ece.cmu.edu/~pueschel/papers/dttalgo.pdf) |
| [R. Gluth, " Regular FFT-related transform kernels for DCT/DST- based polyphase filter banks," Proc. IEEE ICASSP 1991, pp.2205-2208, Toronto, Canada, May 1991.](http://vadkudr.org/Algorithms/REGULAR_FFT_RELATED.zip) |
| [Mathematical relations between discrete trigonometric transform (file in Word-format)](http://vadkudr.org/Algorithms/dtt.zip) |
| [Matlab code presenting some transform relations](http://vadkudr.org/Algorithms/DCT_DST_matlabfiles.zip) |

There are 8 types of discrete cosine and sine transforms. In a paper the first link points Markus Puschel and Jose Moura were developed beautiful theory that combines all 16 trigonometric transforms in one framework and shows that sine and cosine transforms of types 1...4 are related to one another in a very simple fashion. In the same fashion are related to one another sine and cosine transforms of types 5...8. Second link consist paper of Rolf Gluth in that effective approach for computation of sine and cosine transform is presented. Hence, sine and cosine transforms of types 1...4 and Fourier transform (9 most useful transforms) are related to one another in strong mathematical relations. It means that if we have in our DSP-library the code for only one transform (offen only FFT code is available), we can compute any other of 8 transforms using existing code by the small cost of additional pre- and postcalculation. Exhaustive information about this topic you can get visiting [homepage of Markus Puschel](http://www.ece.cmu.edu/~pueschel/publications.html).

Relations between Discrete Trigonometric Transforms\
Matlab reports

|

**DCTI**\
[DCTI<->DSTII](http://vadkudr.org/Algorithms/DTT/DCTI_DSTII/DCTI_DSTII.html)\
[DCTI<->DCTII](http://vadkudr.org/Algorithms/DTT/DCTI_DCTII/DCTI_DCTII.html)\
[DCTI<->DSTI](http://vadkudr.org/Algorithms/DTT/DCTI_DSTI/DCTI_DSTI.html)

 |

**DCTIII**\
[DCTIII<->DCTIV](http://vadkudr.org/Algorithms/DTT/DCTIV_DCTIII/DCTIV_DCTIII.html)\
[DCTIII<->DCTIV(alt)](http://vadkudr.org/Algorithms/DTT/DCTIV_DCTIIIa/DCTIV_DCTIIIa.html)

 |

**DCTV**\
[DCTV<->DSTV](http://vadkudr.org/Algorithms/DTT/DCTV_DSTV/DCTV_DSTV.html)\
[DCTV<->DCTVI](http://vadkudr.org/Algorithms/DTT/DCTV_DCTVI/DCTV_DCTVI.html)\
[DCTV<->DSTVI](http://vadkudr.org/Algorithms/DTT/DCTV_DSTVI/DCTV_DSTVI.html)

 |

**DCTVII**\
[DCTVII<->DSTVII](http://vadkudr.org/Algorithms/DTT/DCTVII_DSTVII/DCTVII_DSTVII.html)\
[DCTVII<->DCTVIII](http://vadkudr.org/Algorithms/DTT/DCTVII_DCTVIII/DCTVII_DCTVIII.html)\
[DCTVII<->DSTVIII](http://vadkudr.org/Algorithms/DTT/DCTVII_DSTVIII/DCTVII_DSTVIII.html)

 |
|

**DSTIII**\
[DSTIII<->DSTIV](http://vadkudr.org/Algorithms/DTT/DSTIV_DSTIII/DSTIV_DSTIII.html)\
[DSTIII<->DSTIV(alt)](http://vadkudr.org/Algorithms/DTT/DSTIV_DSTIIIa/DSTIV_DSTIIIa.html)

 |

**DSTI**\
[DSTI<->DSTII](http://vadkudr.org/Algorithms/DTT/DSTI_DSTII/DSTI_DSTII.html)\
[DSTI<->DCTII](http://vadkudr.org/Algorithms/DTT/DCTII_DSTI/DCTII_DSTI.html)\
[DSTI<->DCTI](http://vadkudr.org/Algorithms/DTT/DCTI_DSTI/DCTI_DSTI.html)

 |

**DSTVII**\
[DCTVII<->DSTVII](http://vadkudr.org/Algorithms/DTT/DCTVII_DSTVII/DCTVII_DSTVII.html)\
[DSTVII<->DCTVIII](http://vadkudr.org/Algorithms/DTT/DCTVIII_DSTVII/DCTVIII_DSTVII.html)\
[DSTVII<->DSTVIII](http://vadkudr.org/Algorithms/DTT/DSTVIII_DSTVII/DSTVIII_DSTVII.html)

 |

**DSTV**\
[DCTV<->DSTV](http://vadkudr.org/Algorithms/DTT/DCTV_DSTV/DCTV_DSTV.html)\
[DSTV<->DCTVI](http://vadkudr.org/Algorithms/DTT/DCTVI_DSTV/DCTVI_DSTV.html)\
[DSTV<->DSTVI](http://vadkudr.org/Algorithms/DTT/DSTVI_DSTV/DSTVI_DSTV.html)

 |
|

**DCTVI**\
[DCTV<->DCTVI](http://vadkudr.org/Algorithms/DTT/DCTV_DCTVI/DCTV_DCTVI.html)\
[DSTV<->DCTVI](http://vadkudr.org/Algorithms/DTT/DCTVI_DSTV/DCTVI_DSTV.html)

 |

**DCTVIII**\
[DCTVII<->DCTVIII](http://vadkudr.org/Algorithms/DTT/DCTVII_DCTVIII/DCTVII_DCTVIII.html)[\
DSTVII<->DCTVIII](http://vadkudr.org/Algorithms/DTT/DCTVIII_DSTVII/DCTVIII_DSTVII.html)

 |

**DCTII**\
[DCTII<->DCTIV](http://vadkudr.org/Algorithms/DTT/DCTIV_DCTII/DCTIV_DCTII.html)\
[DCTII<->DCTIV(alt)](http://vadkudr.org/Algorithms/DTT/DCTIV_DCTIIa/DCTIV_DCTIIa.html)\
[DCTII<->DSTII](http://vadkudr.org/Algorithms/DTT/DCTII_DSTII/DCTII_DSTII.html)\
[DCTII<->DCTI](http://vadkudr.org/Algorithms/DTT/DCTI_DCTII/DCTI_DCTII.html)\
[DCTII<->DSTI](http://vadkudr.org/Algorithms/DTT/DCTII_DSTI/DCTII_DSTI.html)

 |

**DCTIV**\
[DCTIV<->DCTIII](http://vadkudr.org/Algorithms/DTT/DCTIV_DCTIII/DCTIV_DCTIII.html)\
[DCTIV<->DCTIII(alt)](http://vadkudr.org/Algorithms/DTT/DCTIV_DCTIIIa/DCTIV_DCTIIIa.html)\
[DCTIV<->DCTII](http://vadkudr.org/Algorithms/DTT/DCTIV_DCTII/DCTIV_DCTII.html)\
[DCTIV<->DCTII(alt)](http://vadkudr.org/Algorithms/DTT/DCTIV_DCTIIa/DCTIV_DCTIIa.html)

 |
|

**DSTVIII**\
[DCTVII<->DSTVIII](http://vadkudr.org/Algorithms/DTT/DCTVII_DSTVIII/DCTVII_DSTVIII.html)[\
DSTVII<->DSTVIII](http://vadkudr.org/Algorithms/DTT/DSTVIII_DSTVII/DSTVIII_DSTVII.html)

 |

**DSTVI**\
[DCTV<->DSTVI](http://vadkudr.org/Algorithms/DTT/DCTV_DSTVI/DCTV_DSTVI.html)[\
DSTV<->DSTVI](http://vadkudr.org/Algorithms/DTT/DSTVI_DSTV/DSTVI_DSTV.html)

 |

**DSTIV**\
[DSTIV<->DSTIII](http://vadkudr.org/Algorithms/DTT/DSTIV_DSTIII/DSTIV_DSTIII.html)\
[DSTIV<->DSTIII(alt)](http://vadkudr.org/Algorithms/DTT/DSTIV_DSTIIIa/DSTIV_DSTIIIa.html)\
[DSTIV<->DSTII](http://vadkudr.org/Algorithms/DTT/DSTIV_DSTII/DSTIV_DSTII.html)\
[DSTIV<->DSTII(alt)](http://vadkudr.org/Algorithms/DTT/DSTIV_DSTIIa/DSTIV_DSTIIa.html)

 |

**DSTII**\
[DSTII<->DSTIV](http://vadkudr.org/Algorithms/DTT/DSTIV_DSTII/DSTIV_DSTII.html)\
[DSTII<->DSTIV(alt)](http://vadkudr.org/Algorithms/DTT/DSTIV_DSTIIa/DSTIV_DSTIIa.html)\
[DSTII<->DCTII](http://vadkudr.org/Algorithms/DTT/DCTII_DSTII/DCTII_DSTII.html)\
[DSTII<->DCTI](http://vadkudr.org/Algorithms/DTT/DCTI_DSTII/DCTI_DSTII.html)\
[DSTII<->DSTI](http://vadkudr.org/Algorithms/DTT/DSTI_DSTII/DSTI_DSTII.html)

 |

Additional notes:

-   May be it is interesting for someone - I placed here matlab scripts [GenDTT.zip](http://vadkudr.org/Algorithms/DTT/GenDTT.zip) for generating transform generations. Just place proper pair of transforms at the line 143 and run. Just be sure that transform pair is from cells of the same color.
-   I made detailed comments to the DST2<->DCT1 transform conversion example presented in the M. Pueschel's paper mentioned above ([DTT transform relationship derivation(detailed).pdf](http://vadkudr.org/Algorithms/DTT/DTT%20transform%20relationship%20derivation%28detailed%29.pdf)).
