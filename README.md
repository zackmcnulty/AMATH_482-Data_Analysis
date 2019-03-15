# AMATH 482: Data Analysis and Signal Processing
- video lectures available at : https://faculty.washington.edu/kutz/KutzBook/KutzBook.html
- If you develop a new technique, show it works on a small system where you KNOW the truth.
- Machine learning is cool and all, but it should not be the first tool you go for. There
are a lot simpler methods (SVD, Fourier Transform, etc) that might do the trick or make
the problem easier.
- Solving a PDE? Convert it to an ODE and use ode45. Consider taking the Fourier transform to convert a spatial temporal
 signal to eliminate the spatial term.

# Lecture Topics

### Data Analysis
- Fourier Transforms
- Averaging and filtering in frequency domain
- Short-term transforms (Gabor Transforms)
- Wavelet Transforms (no smooth filters)

### Image Processing
- Singular Value Decomposition (SVD)
- Principal Component Analysis (PCA)
- Low rank (dimension) approximations to data
- Dimensionality analysis of data


### Classification & Clustering
##### Unsupervised Learning
- K-means
- Hierarchial Clustering
- Gaussian Mixture Models

##### Supervised Learning
- K-nearest Neighbors
- Linear Discriminant Analysis (pg 173 of pdf)
- Support Vector Machines (pg 179 of pdf)
- Classification & Regression Trees (pg 184 of pdf)


### Dynamic Mode Decomposition
- http://dmdbook.com/
- Nonorthogonal Mode decomposition
- Extracting spatio-temporal properties


### Compressive Sensing
- Sensing in large systems when you only have sparse measurements

### Reduced Order Modeling
- Solving PDEs by reducing to ODEs 
- Benefits of Fourier Transform in PDE --> ODE conversion


# Homework Topics

### HW 1
Using ultrasound data to track the location of a foreign object throughout
a noisy medium. Techniques for handeling zero-mean frequency-based noise.

* Fourier Transforms
* Averaging in the Frequency Domain
* Filtering for Frequencies


### HW 2
- Short-term fourier transforms (Gabor transforms)
- Spectrograms
- Music decomposition

### HW 3
- Singular Value Decomposition (SVD)
- Dimensionality Analysis and Low-Rank (dimension) approximations of data
- Principal Component Analysis (PCA)

### HW 4
- Data Preprocessing
- Classification and Clustering

### HW 5

- Dynamic Mode decomposition
- Separating image foreground from background
- Extracting moving objects from video
