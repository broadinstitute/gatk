import $ivy.`ch.qos.logback:logback-core:1.4.5`
import $ivy.`ch.qos.logback:logback-classic:1.4.5`
import $ivy.`com.azure:azure-cosmos:4.41.0`
import $ivy.`com.lihaoyi:ammonite-ops_2.13:2.4.1`
import $ivy.`org.slf4j:slf4j-api:2.0.6`
import $ivy.`io.projectreactor.tools:blockhound:1.0.7.RELEASE`


import com.azure.cosmos.CosmosAsyncContainer
import com.azure.cosmos.CosmosException
import com.azure.cosmos.implementation.HttpConstants
import com.azure.cosmos.models.CosmosBulkExecutionOptions
import com.azure.cosmos.models.CosmosBulkItemResponse
import com.azure.cosmos.models.CosmosBulkOperationResponse
import com.azure.cosmos.models.CosmosItemOperation
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import reactor.core.publisher._
import reactor.core.scheduler.Schedulers
import java.util.concurrent.Semaphore


object BulkWriter {
  private val logger = LoggerFactory.getLogger(classOf[BulkWriter])
}

class BulkWriter(private val cosmosAsyncContainer: CosmosAsyncContainer) {
  final private val bulkInputEmitter: Sinks.Many[CosmosItemOperation] = Sinks.many.unicast.onBackpressureBuffer
  final private val cpuCount = Runtime.getRuntime.availableProcessors
  BulkWriter.logger.info(f"${cpuCount} CPUs available")
  //Max items to be buffered to avoid out of memory error
  final private val maxItems = 1024 * 167 / cpuCount
  BulkWriter.logger.info(f"Buffering for ${maxItems} items max")
  final private val semaphore = new Semaphore(maxItems)
  final private val emitFailureHandler: Sinks.EmitFailureHandler = (signalType: SignalType, emitResult: Sinks.EmitResult) => {
    if (emitResult == Sinks.EmitResult.FAIL_NON_SERIALIZED) {
      BulkWriter.logger.info("emitFailureHandler - Signal: [{}], Result: [{}]", signalType, emitResult)
      true
    }
    else {
      BulkWriter.logger.error("emitFailureHandler - Signal: [{}], Result: [{}]", signalType, emitResult)
      false
    }
  }

  def scheduleWrites(cosmosItemOperation: CosmosItemOperation): Unit = {
    BulkWriter.logger.info("in scheduleWrites...")
    while (!semaphore.tryAcquire) BulkWriter.logger.info("Unable to acquire permit")
    BulkWriter.logger.info("Acquired permit")
    scheduleInternalWrites(cosmosItemOperation)
  }

  private def scheduleInternalWrites(cosmosItemOperation: CosmosItemOperation): Unit = {
    bulkInputEmitter.emitNext(cosmosItemOperation, emitFailureHandler)
  }

  def execute: Flux[CosmosBulkOperationResponse[_]] = this.execute(null)

  def execute(bulkOptions: CosmosBulkExecutionOptions): Flux[CosmosBulkOperationResponse[_]] = {
    cosmosAsyncContainer.executeBulkOperations(bulkInputEmitter.asFlux, if (bulkOptions == null) new CosmosBulkExecutionOptions else bulkOptions).publishOn(Schedulers.boundedElastic).map((bulkOperationResponse: CosmosBulkOperationResponse[AnyRef]) => {
      processBulkOperationResponse(bulkOperationResponse.getResponse, bulkOperationResponse.getOperation, bulkOperationResponse.getException)
      semaphore.release()
      bulkOperationResponse

    })
  }

  private def processBulkOperationResponse(itemResponse: CosmosBulkItemResponse, itemOperation: CosmosItemOperation, exception: Exception): Unit = {
    if (exception != null) handleException(itemOperation, exception)
    else processResponseCode(itemResponse, itemOperation)
  }

  private def processResponseCode(itemResponse: CosmosBulkItemResponse, itemOperation: CosmosItemOperation): Unit = {
    if (itemResponse.isSuccessStatusCode) BulkWriter.logger.info("The operation for Item ID: [{}]  Item PartitionKey Value: [{}] completed successfully " + "with a response status code: [{}]", itemOperation.getId, itemOperation.getPartitionKeyValue, itemResponse.getStatusCode)
    else if (shouldRetry(itemResponse.getStatusCode)) {
      BulkWriter.logger.info("The operation for Item ID: [{}]  Item PartitionKey Value: [{}] will be retried", itemOperation.getId, itemOperation.getPartitionKeyValue)
      //re-scheduling
      scheduleWrites(itemOperation)
    }
    else BulkWriter.logger.info("The operation for Item ID: [{}]  Item PartitionKey Value: [{}] did not complete successfully " + "with a response status code: [{}]", itemOperation.getId, itemOperation.getPartitionKeyValue, itemResponse.getStatusCode)
  }

  private def handleException(itemOperation: CosmosItemOperation, exception: Exception): Unit = {
    if (!exception.isInstanceOf[CosmosException]) BulkWriter.logger.info("The operation for Item ID: [{}]  Item PartitionKey Value: [{}] encountered an unexpected failure", itemOperation.getId, itemOperation.getPartitionKeyValue)
    else if (shouldRetry(exception.asInstanceOf[CosmosException].getStatusCode)) {
      BulkWriter.logger.info("The operation for Item ID: [{}]  Item PartitionKey Value: [{}] will be retried", itemOperation.getId, itemOperation.getPartitionKeyValue)
      //re-scheduling
      scheduleWrites(itemOperation)
    }
  }

  private def shouldRetry(statusCode: Int) = statusCode == HttpConstants.StatusCodes.REQUEST_TIMEOUT || statusCode == HttpConstants.StatusCodes.TOO_MANY_REQUESTS
}
