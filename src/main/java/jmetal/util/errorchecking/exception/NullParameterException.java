package jmetal.util.errorchecking.exception;

@SuppressWarnings("serial")
public class NullParameterException extends RuntimeException {
  public NullParameterException() {
    super("The parameter is null") ;
  }
}
