import org.ojalgo.OjAlgoUtils;
import org.ojalgo.RecoverableCondition;
import org.ojalgo.matrix.PrimitiveMatrix;
import org.ojalgo.matrix.decomposition.QR;
import org.ojalgo.matrix.store.ElementsSupplier;
import org.ojalgo.matrix.store.MatrixStore;
import org.ojalgo.matrix.store.PhysicalStore;
import org.ojalgo.matrix.store.PrimitiveDenseStore;
import org.ojalgo.matrix.task.InverterTask;
import org.ojalgo.matrix.task.SolverTask;
import org.ojalgo.netio.BasicLogger;
import org.ojalgo.random.Normal;

/**
 * Getting Started with ojAlgo
 */
public class GettingStarted {

	public static void main(String[] args) {

		BasicLogger.debug();
		BasicLogger.debug(GettingStarted.class);
		BasicLogger.debug(OjAlgoUtils.getTitle());
		BasicLogger.debug(OjAlgoUtils.getDate());
		BasicLogger.debug();

		PrimitiveMatrix.Factory matrixFactory = PrimitiveMatrix.FACTORY;
		PhysicalStore.Factory<Double, PrimitiveDenseStore> storeFactory = PrimitiveDenseStore.FACTORY;
		// PrimitiveMatrix.Factory and PrimitiveDenseStore.Factory are very similar.
		// Every factory in ojAlgo that makes 2D-structures
		// extends/implements the same interface.

		PrimitiveMatrix matrixA = matrixFactory.makeEye(5, 5);
		// Internally this creates an "eye-structure" - not a large array.
		PrimitiveDenseStore storeA = storeFactory.makeEye(5, 5);
		// A PrimitiveDenseStore is always a "full array". No smart data
		// structures here.

		PrimitiveMatrix matrixB = matrixFactory.makeFilled(5, 3, new Normal(5.0, 2.0));
		PrimitiveDenseStore storeB = storeFactory.makeFilled(5, 3, new Normal(5.0, 2.0));
		// When you create a matrix with random elements you can specify
		// their distribution.

		/* Matrix multiplication */

		PrimitiveMatrix matrixC = matrixA.multiply(matrixB);
		// Multiplying two PrimitiveMatrix:s is trivial. There are no
		// alternatives, and the returned product is a PrimitiveMatrix
		// (same as the inputs).

		// Doing the same thing using PrimitiveDenseStore (MatrixStore) you
		// have options...

		BasicLogger.debug("Different ways to do matrix multiplication with " + "MatrixStore:s");
		BasicLogger.debug();

		MatrixStore<Double> storeC = storeA.multiply(storeB);
		// One option is to do exactly what you did with PrimitiveMatrix.
		// The only difference is that the return type is MatrixStore rather
		// than PhysicalStore, PrimitiveDenseStore or whatever else you input.
		BasicLogger.debug("MatrixStore MatrixStore#multiply(MatrixStore)", storeC);

		PrimitiveDenseStore storeCPreallocated = storeFactory.makeZero(5, 3);
		// Another option is to first create the matrix that should hold the
		// resulting product,
		storeA.multiply(storeB, storeCPreallocated);
		// and then perform the multiplication. This enables reusing memory
		// (the product matrix).
		BasicLogger.debug("void MatrixStore#multiply(Access1D, ElementsConsumer)", storeCPreallocated);

		ElementsSupplier<Double> storeCSupplier = storeB.premultiply(storeA);
		// A third option is the premultiply method:
		// 1) The left and right argument matrices are interchanged.
		// 2) The return type is an ElementsSupplier rather than
		//    a MatrixStore.
		// This is because the multiplication is not yet performed.
		// It is possible to define additional operation on
		// an ElementsSupplier.
		MatrixStore<Double> storeCLater = storeCSupplier.get();
		// The multiplication, and whatever additional operations you defined,
		// is performed when you call #get().
		BasicLogger.debug("ElementsSupplier MatrixStore#premultiply(Access1D)", storeCLater);

		// A couple of more alternatives that will do the same thing.
		storeCPreallocated.fillByMultiplying(storeA, storeB);
		BasicLogger.debug("void ElementsConsumer#fillByMultiplying(Access1D, Access1D)", storeCLater);
		storeCSupplier.supplyTo(storeCPreallocated);
		BasicLogger.debug("void ElementsSupplier#supplyTo(ElementsConsumer)", storeCLater);

		/* Inverting matrices */

		matrixA.invert();
		// If you want to invert a PrimitiveMatrix then just do it...
		// The matrix doesn't even have to square or anything.
		// You'll get a pseudoinverse or whatever is possible.

		// With MatrixStore:s you need to use an InverterTask
		InverterTask<Double> inverter = InverterTask.PRIMITIVE.make(storeA);
		// There are many implementations of that interface. This factory
		// method will return one that may be suitable, but most likely you
		// will want to choose implementation based on what you know about
		// the matrix.
		try {
			inverter.invert(storeA);
		} catch (RecoverableCondition e) {
			// Will throw and exception if inversion fails, rethrowing it.
			throw new RuntimeException(e);
		}

		/* Solving equation system */

		matrixA.solve(matrixC);
		// PrimitiveMatrix "is" the equation system body.
		// You just supply the RHS, and you get the solution.
		// If necessary it will be a least squares solution, or something
		// derived from the pseudoinverse.

		SolverTask<Double> solver = SolverTask.PRIMITIVE.make(storeA, storeC);
		// Similar to InverterTask:s there are SolverTask:s for the lower level stuff
		try {
			solver.solve(storeA, storeC);
		} catch (RecoverableCondition e) {
			// Will throw and exception if solving fails, rethrowing it.
			throw new RuntimeException(e);
		}

		// Most likely you want to do is to instantiate some matrix
		// decomposition (there are several).

		QR<Double> qr = QR.PRIMITIVE.make(storeA);
		// You supply a typical matrix to the factory to allow it to choose implementation
		// (depending on the size/shape).
		qr.decompose(storeA);
		// Then you do the decomposition
		if (qr.isSolvable()) {
			// You should verify that the equation system is solvable,
			// and do something else if it is not.
			qr.getSolution(storeC);
		} else {
			throw new RuntimeException("Cannot solve the equation system");
		}

		/* Setting individual elements */

		storeA.set(3, 1, 3.14);
		storeA.set(3, 0, 2.18);
		// PhysicalStore instances are naturally mutable. If you want to set
		// or modify something - just do it

		PrimitiveMatrix.DenseReceiver matrixBuilder = matrixA.copy();
		// PrimitiveMatrix is immutable. To modify anything, you need to copy
		// it to builder instance.

		matrixBuilder.add(3, 1, 3.14);
		matrixBuilder.add(3, 0, 2.18);

		matrixA = matrixBuilder.build();

		/* Creating matrices by explicitly setting all elements */

		double[][] data = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};

		matrixFactory.rows(data);
		storeFactory.rows(data);

		// If you don't want/need to first create some (intermediate) array
		// for the elements, you can of course set them on the matrix
		// directly.
		PrimitiveDenseStore storeZ = storeFactory.makeEye(3, 3);

		// Since PrimitiveMatrix is immutable this has to be done via
		// a builder.
		PrimitiveMatrix.DenseReceiver matrixZBuilder = matrixFactory.makeDense(3, 3);

		for (int j = 0; j < 3; j++) {
			for (int i = 0; i < 3; i++) {
				matrixZBuilder.set(i, j, i * j);
				storeZ.set(i, j, i * j);
			}
		}

		PrimitiveMatrix matrixZ = matrixZBuilder.get();

		BasicLogger.debug();
		BasicLogger.debug("PrimitiveMatrix Z", matrixZ);
	}
}