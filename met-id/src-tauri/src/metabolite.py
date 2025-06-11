import sys
import grpc
import asyncio
import streaming_pb2
import streaming_pb2_grpc
from rdkit import Chem

def substruct_matches(smiles_strings, smarts_pattern):
    smarts_mol = Chem.MolFromSmarts(smarts_pattern)
    if smarts_mol is None:
        return "Invalid SMARTS pattern"
    
    results = []
    for smiles in smiles_strings:
        input_mol = Chem.MolFromSmiles(smiles)
        if input_mol is None:
            results.append(0)
        else:
            matches = len(input_mol.GetSubstructMatches(smarts_mol))
            results.append(matches)
    #print("results:", results)
    return results


async def perform_computation(smiles, smarts):
    """Simulates a time-consuming computation."""
    #print(f"Starting computation for: {message}")
    #print("smarts: ", smarts, "smiles:", smiles)
    result = substruct_matches(smiles, smarts)
    #await asyncio.sleep(0.1)  # Simulate computation delay

    #print(f"Finished computation for: {result}")
    return ", ".join(str(x) for x in result)

async def stream_batched_messages(stub, smarts):
    
    async def request_stream():
        rcvd = 1
        done = 1

        while True:
            # Wait for the server's response
            #print("Waiting for server response...")
            received = await responses_queue.get()
            if received is None:
                #print("No more messages to process. Ending stream.")
                break

            #print(f"Received from server: {received.messages}")

            # If server sends "stop", send a stop response and break
            if "stop" in received.messages:
                #print("Received stop signal from server. Sending stop response...")
                yield streaming_pb2.MessageBatch(messages=["stop"])
                break

            # Perform computation and send back the result
            #for message in received.messages:
            print("received: ", rcvd)
            rcvd += 1
            result = await perform_computation(received.messages, smarts)  # Compute the result
            print("sent: ", done)
            done += 1
            #print(f"Sending result back to server: {result}")
            yield streaming_pb2.MessageBatch(messages=[result])

    responses_queue = asyncio.Queue()

    #print("Starting bidirectional streaming...")
    try:
        async for response in stub.StreamBatchedMessages(request_stream()):
            await responses_queue.put(response)
    except grpc.aio.AioRpcError as e:
        print(f"gRPC error: {e}")
    finally:
        await responses_queue.put(None)
        #print("Stream completed.")

async def main(smarts):
    async with grpc.aio.insecure_channel('localhost:50051') as channel:
        stub = streaming_pb2_grpc.StreamingServiceStub(channel)
        #print("Client started.")
        await stream_batched_messages(stub, smarts)

if __name__ == '__main__':
    smarts = sys.argv[1]
    asyncio.run(main(smarts))
